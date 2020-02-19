import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import spearmanr, linregress, pearsonr, rankdata
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def import_survival_data(filepath):
    """ Import survival data either having two columns (time, survival) or one
    column (time)

    Args:
        filepath (str): path to survival data file

    Returns:
        pd.DataFrame: returned data frame
    """
    with open(filepath, 'r') as f:
        cols = len(f.readline().split(','))
    if cols == 2:
        df = pd.read_csv(filepath, sep=',', header=0, names=['Time', 'Survival'])
    elif cols == 1:
        #TODO remove repeating points at the end of the tail
        df = pd.read_csv(filepath, sep=',', header=None, names=['Time'])
        df.loc[:, 'Survival'] = np.linspace(0, 100, num=df.shape[0])

    df = df.drop_duplicates()
    return df


def parse_input(filepath):
    """ Parse input file containing names and file paths

    Args:
        filepath (str): path to input file

    Returns:
        list: list containing name and file paths
    """
    ls = []
    with open(filepath) as f:
        for line in f.readlines():
            tokens = line.strip().split(',')
            ls.append((tokens[0], import_survival_data(tokens[1])))
    return ls


def add_pseudo_point(ori_df):
    """ Add pseudo point (end time, survival 0.01) for scaling purposes

    Args:
        ori_df (pd.DataFrame): survival data

    Returns:
        pd.DataFrame: copy of dataframe with pseudo point added
    """
    # add pseudo point with survival (end time, 0) for scaling purposes
    df = ori_df.copy()
    df.index = df.index + 1
    df.loc[0, :] = [df.at[1, 'Time'], 0.01]
    df = df.sort_index()
    return df


def remove_pseudo_point(ori_df):
    """ Remove pseudo point (end time, survival 0.01) added for scaling purposes

    Args:
        ori_df (pd.DataFrame): survival data

    Returns:
        pd.DataFrame: copy of dataframe with pseudo point removed
    """
    df = ori_df.copy()
    df = df.drop(index=0)
    df.index = df.index - 1
    return df


def interpolate(df, x='Time', y='Survival', kind='linear'):
    return interp1d(df[x], df[y], kind=kind, fill_value='extrapolate')


def _comb_prob_theor(p_high, p_low, rho):
    """ Theoretical approcach. When P(A) > P(B), P(A+B) = P(A) + (1-P(A)) * P(B) * (1-rho)

    Args:
        p_high (np.ndarray): 1-D array of higher probability of survival at each timepoint
        p_low (np.ndarray): 1-D array of lower probability of survival at each timepoint
        rho (float): correlation coefficient

    Returns:
        np.ndarray: 1-D array of predicted survival (%) at each timepoint for
                    combination therapy
    """
    return (p_high + ((1 - p_high) * p_low * (1 - rho))) * 100


def combined_prob(df_a, df_b, timecourse, rho=0.3):
    """ Calculate predicted survival during time course for combination therapy
    based on theoretical approach.

    Args:
        df_a (pd.DataFrame): survival data for treatment A
        df_b (pd.DataFrame): survival data for treatment B
        timecourse (list-like obj): timepoints to predict survival
        rho (float): correlation coefficient (default: 0.3)

    Returns:
        np.ndarray: 1-D array of predicted survival (%) at each timepoint for
                    combination therapy
    """
    # interpolate
    f_a = interpolate(df_a, x='Time', y='Survival')
    f_b = interpolate(df_b, x='Time', y='Survival')
    # convert percentage to decimals
    prob_a, prob_b = f_a(timecourse) / 100, f_b(timecourse) / 100
    return _comb_prob_theor(np.fmax(prob_a, prob_b), np.fmin(prob_a, prob_b), rho)


def sample_indep_response(df, patients, n=2000):
    """ Randomly sample PFS response (time) for n-patients.

    Args:
        df (pd.DataFrame): survival data
        n (int): number of points (patients) (default: 2000)
        patients (np.ndarray): 1-D array of n-patients from in 0-100 scale.

    Returns:
        np.ndarray: 1-D array of randomly sampled PFS time for n-patients
    """
    rand_idx = np.random.permutation(n)
    f = interpolate(df, x='Survival', y='Time')
    survive_time = f(patients)
    return survive_time[rand_idx]


def fit_rho(n, desired_rho, size=200):
    """ Shuffle index of one data to make two datasets to have desired spearmen
    correlation coefficient (random sampling version)
    Args:
        n (int): number of data points
        size (int): generate random integer [0, size) to shuffle index
        desired_rho (float): desired spearman correlation coefficient

    Returns:
        np.ndarray: shuffled index
    """
    if size < 0:
        size += 10

    permit = 0.005
    shuffled = [i + np.random.randint(size) for i in range(n)]
    rho, _ = spearmanr(range(n), shuffled)
    if rho < desired_rho - permit:
        return fit_rho(n, desired_rho, size - 5)
    elif rho > desired_rho + permit:
        return fit_rho(n, desired_rho, size + 5)
    else:
        return np.argsort(shuffled)


def fit_rho2(n, rho):
    """ Shuffle indices of two sorted dataset to make two datasets to have desired Spearman
    correlation coefficient (copulas version)
    Based on: https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables/313138#313138

    Args:
        n (int): number of data points
        rho (float): desired spearman correlation coefficient

    Returns:
        tuple: tuple of shuffled indices (np.ndarray)
    """
    rho_p = 2 * np.sin(rho * np.pi / 6)  # convert to pearson's corr coef
    x = np.random.permutation(n)
    y = np.random.permutation(n)
    result = linregress(x, y)
    residuals = y - (result[0] * x + result[1])
    new_y = rho_p * np.std(residuals) * x + np.sqrt(1 - rho_p**2) * np.std(x) * residuals
    return ((rankdata(x) - 1).astype(int), (rankdata(new_y) - 1).astype(int))


def sample_joint_response(ori_a, ori_b, patients, rho=0, n=2000):
    """ Calculate predicted PFS time for n-patients in combination therapy based
    on sampling approach.

    Args:
        ori_a (pd.DataFrame): survival data for treatment A
        ori_b (pd.DataFrame): survival data for treatment B
        patients (np.ndarray): 1-D array of n-patients from in 0-100 scale.
        rho (float): desired spearman correlation coefficient (default: 0)
        n (int): number of data points (patients) (default: 2000)

    Returns:
        list: list of predicted PFS for n-patients in combination therapy
    """
    df_a, df_b = ori_a, ori_b
    if rho == 0:
        return sorted(np.maximum(sample_indep_response(df_a, patients, n=n),
                                 sample_indep_response(df_b, patients, n=n)), reverse=True)
    else:
        fa = interpolate(df_a, x='Survival', y='Time')
        fb = interpolate(df_b, x='Survival', y='Time')
        #shuffled = fit_rho2(n, rho)
        idx_a, idx_b = fit_rho2(n, rho)
        return sorted(np.maximum(fa(patients)[idx_a], fb(patients)[idx_b]), reverse=True)


def adjust_response(df, time, response):
    """ Adjust response to make survival at {time} to {response} survival.

    Args:
        df (pd.DataFrame): survival data frame
        time (float):
        response (float):

    Returns:
        pd.DataFrame: updated survival data frame
    """
    df = add_pseudo_point(df)
    cutoff_idx = df[df['Time'] < time]['Survival'].idxmin()
    res = df.loc[:cutoff_idx, :].copy()
    nonres = df.loc[cutoff_idx:, :].copy()
    res.loc[:, 'Survival'] = np.interp(res['Survival'],
                                       (res['Survival'].min(),
                                        res['Survival'].max()),
                                       (0, response))
    nonres.loc[:, 'Survival'] = np.interp(nonres['Survival'],
                                          (nonres['Survival'].min(),
                                           nonres['Survival'].max()),
                                          (response, 100))
    df = res.append(nonres)
    df = remove_pseudo_point(df)
    return df


def median_pfs(df_list, ax):
    """ Plots dashed lines for median PFS.
    Args:
        df_list (list): list of survival data frames
        ax (matplotlib.axes.Axes): axis to plot
    """
    med_pfs_list = []
    # calculate median PFS
    for df in df_list:
        inv_f = interpolate(df, x='Survival', y='Time')
        med_pfs_list.append(inv_f(50))

    ax.hlines(50, 0, max(med_pfs_list), linestyle='--', linewidth=1)

    for med_pfs in med_pfs_list:
         ax.vlines(med_pfs, 50, 0, linestyle='--', linewidth=1)


def main():
    parser = argparse.ArgumentParser(
        description="Predict combination effect of drugs",
        epilog='Implemented by Haeun Hwangbo')
    parser.add_argument(
        'input', help='CSV File describing name and path of the survival data.The ordering of the paths should be drug A, B, A+B. (Column 1: name of drug, Column 2: path to survival data)')
    parser.add_argument('--min-rho', type=float,
                        default=0.1, help='Minimum rho (default: 0.1)')
    parser.add_argument('--max-rho', type=float,
                        default=0.5, help='Maximum rho (default: 0.5)')
    parser.add_argument('--adj-respA', nargs=2, type=float,
                        help='Adjust resposne of treatment A to (survival) at (time). Should be in the format: --adj-respA time survival')
    parser.add_argument('--adj-respB', nargs=2, type=float,
                        help='Adjust resposne of treatment B to (survival) at (time). Should be in the format: --adj-respB time survival')
    parser.add_argument('--fig-height', type=int, default=5,
                        help='Output figure height (default: 5)')
    parser.add_argument('--fig-width', type=int, default=10,
                        help='Output figure width (default: 10)')
    parser.add_argument(
        '--out-prefix', help='Output file prefix. If not specified, it will be the names of the drugs.')
    parser.add_argument('--extension', default='pdf', choices=['pdf', 'png', 'jpg'],
                        help='File extension for output figure (default: pdf)')
    parser.add_argument('--predict-type', default='theor', choices=['theor', 'stoch'],
                        help='How to predict combination effect. Calculate by equation or stochastic sampling (default: theor)')
    parser.add_argument('--N', type=int, default=5000,
                        help='Number of data points to use in prediction (default: 5000)')
    args = parser.parse_args()

    # get input
    treatments = parse_input(args.input)

    name_a, df_a = treatments[0]
    name_b, df_b = treatments[1]
    name_ab, df_ab = treatments[2]

    if args.adj_respA:
        df_a = adjust_response(df_a, args.adj_respA[0], args.adj_respA[1])
    if args.adj_respB:
        df_b = adjust_response(df_b, args.adj_respB[0], args.adj_respB[1])

    # predict combination effect
    N = args.N
    if args.predict_type == 'theor':
        timepoints = np.linspace(
            0, min(df_a['Time'].max(), df_b['Time'].max()), num=N)

        predicted = pd.DataFrame({"Time": timepoints,
                                  "Survival": combined_prob(df_a, df_b, timepoints,
                                                            rho=(args.min_rho + args.max_rho) / 2)})
    else:  # stoch
        patients = np.linspace(0, 100, num=N)
        predicted = pd.DataFrame({'Survival': patients,
                                  'Time': sample_joint_response(df_a, df_b, patients, n=N,
                                                                rho=(args.min_rho + args.max_rho) / 2)})
    # plot survival curve
    fig, ax=plt.subplots(figsize = (args.fig_width, args.fig_height))
    sns.despine()
    sns.lineplot(x = 'Time', y = 'Survival', data = df_a,
                 label = name_a, ax = ax)
    sns.lineplot(x = 'Time', y = 'Survival', data = df_b,
                 label = name_b, ax = ax)
    sns.lineplot(x = 'Time', y = 'Survival', data = df_ab,
                 label = name_ab, ax = ax)
    sns.lineplot(x = 'Time', y = 'Survival', data = predicted,
                 color = 'black', label = 'Combination Predicted', ax = ax)
    if args.predict_type == 'theor':
        ax.fill_between(timepoints,
                        combined_prob(df_a, df_b, timepoints,
                                      rho=args.min_rho),
                        combined_prob(df_a, df_b, timepoints,
                                      rho=args.max_rho),
                        alpha = 0.3, color = 'gray')
    else:
        ax.fill_betweenx(patients,
                         sample_joint_response(
                             df_a, df_b, patients, n=N, rho=args.min_rho),
                         sample_joint_response(
                             df_a, df_b, patients, n=N, rho=args.max_rho),
                         alpha = 0.3, color = 'gray')
    median_pfs([df_a, df_b, df_ab, predicted], ax)
    ax.set_xlim(0, max([df_a['Time'].max(), df_b['Time'].max(), df_ab['Time'].max()]) + 1)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel('Survival (%)')
    fig.tight_layout()
    # save output figure
    if args.out_prefix is None:
        fig.savefig('./{0}_{1}_combination_kmplot.pdf'.format(name_a, name_b))
    else:
        fig.savefig(args.out_prefix + '.' + args.extension)


if __name__ == '__main__':
    main()
