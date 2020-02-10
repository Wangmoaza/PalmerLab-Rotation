import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def import_survival_data(filepath):
    df = pd.read_csv(filepath, sep=',', header=0, names=['Time', 'Survival'])
    df.loc[:, 'Survival'] = df['Survival']
    df = df.drop_duplicates()
    return df


def _comb_prob(p_high, p_low, rho):
    # When P(A) > P(B)
    # P(A+B) = P(A) + (1-P(A)) * P(B) * rho
    return (p_high + ((1 - p_high) * p_low * (1 - rho))) * 100


def combined_prob(f_a, f_b, timecourse, rho=0.3):
    # convert percentage to decimals
    prob_a, prob_b = f_a(timecourse) / 100, f_b(timecourse) / 100
    return _comb_prob(np.fmax(prob_a, prob_b), np.fmin(prob_a, prob_b), rho)


def interpolate(df, kind='linear'):
    return interp1d(df['Time'], df['Survival'], kind=kind, fill_value='extrapolate')


def parse_input(filepath):
    ls = []
    with open(filepath) as f:
        for line in f.readlines():
            tokens = line.strip().split(',')
            ls.append((tokens[0], import_survival_data(tokens[1])))
    return ls


def main():
    parser = argparse.ArgumentParser(
        description="Predict combination effect of drugs",
        epilog='Implemented by Haeun Hwangbo')
    parser.add_argument(
        'input', help='CSV File describing name and path of the survival data.\nThe ordering of the files should be drug A, B, A+B.\n(Column 1: name of drug, Column 2: path to survival data)')
    parser.add_argument('--min-rho', type=float,
                        default=0.1, help='Minimum rho (default: 0.1)')
    parser.add_argument('--max-rho', type=float,
                        default=0.5, help='Maximum rho (default: 0.5)')
    parser.add_argument('--fig-height', type=int, default=5,
                        help='Output figure height (default: 5)')
    parser.add_argument('--fig-width', type=int, default=10,
                        help='Output figure width (default: 10)')
    parser.add_argument('--out-prefix', help='Output file prefix. If not specified, it will be the names of the drugs.')
    parser.add_argument('--extension', default = 'pdf', choices=['pdf', 'png', 'jpg'],
                        help='File extension for output figure')
    args = parser.parse_args()

    treatments = parse_input(args.input)

    name_a, df_a = treatments[0]
    name_b, df_b = treatments[1]
    name_ab, df_ab = treatments[2]

    # interpolate
    f_a = interpolate(df_a)
    f_b = interpolate(df_b)

    timepoints = np.arange(
        0, int(min(df_a['Time'].max(), df_b['Time'].max())), 0.01)

    predicted = pd.DataFrame({"Time": timepoints,
                              "Survival": combined_prob(f_a, f_b, timepoints,
                                                        rho=(args.min_rho + args.max_rho) / 2)})

    # draw survival curves
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    sns.despine()
    sns.lineplot(x='Time', y='Survival', data=df_a,
                 label=name_a, ax=ax)
    sns.lineplot(x='Time', y='Survival', data=df_b,
                 label=name_b, ax=ax)
    sns.lineplot(x='Time', y='Survival', data=df_ab,
                 label=name_ab, ax=ax)
    sns.lineplot(x='Time', y='Survival', data=predicted,
                 color='black', label='Combination Predicted', ax=ax)
    ax.fill_between(timepoints,
                    combined_prob(f_a, f_b, timepoints, rho=args.min_rho),
                    combined_prob(f_a, f_b, timepoints, rho=args.max_rho),
                    alpha=0.3, color='gray')
    ax.set_xlim(0)
    ax.set_ylim(0)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel('Survival (%)')
    # save output figure
    if args.out_prefix is None:
        fig.savefig('./{0}_{1}_combination_kmplot.pdf'.format(name_a, name_b))
    else:
        fig.savefig(args.out_prefix + '.' + args.extension)


if __name__ == '__main__':
    main()
