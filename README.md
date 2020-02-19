# PalmerLab-Rotation

## Rotation with Dr. Adam Palmer

## Lab Notes

### 2020-02-07

1.  Worked on predicting combination therapy effect
    -   check output with jupyter notebook

### 2020-02-10

1.  Stand-alone program for predicting combination therapy effect.

### 2020-02-12

1.  Added median PFS line to plot function

### 2020-02-17

1.  Implemented sampling version of prediction

### 2020-02-18

1.  Found bug in sampling version
2.  Changed way of generating joint uniform distribution with desired rho
    -   randomly adding integer to the second data was taking too long
    -   use copulas equation to generate joint distribution
    -   faster algorithm

### 2020-02-19

1.  Fixed bug in sampling version

    1.  np.argsort wasn't what I expected it to do (ranking) -> changed it to scipy.stats.rankdata
    2.  need to put patient 0 - 100 not max(min(data1), min(data2)) when generating joint distribution

2.  Copulas equation

    -   copulas is using pearson's rho, but we want to set the spearman rho
        -   so convert the spearman rho to the equivalent pearson's rho

3.  added N parameter

4. Confirmed that equation-based and sampling-based methods match

5. Re-run on adjusted resposne
   - time point: 1.5, 3, 3.5
   - adjusted survival at time point: 60, 70, 80, 90
   - N: 50,000
   - method: equation, sampling
