#!/usr/bin/env bash

#for time in 1.5 3 4.5
for time in 3
do
#    for response in 60 70 80 90
    for response in 58.5
    do
        python predict_comb_therapy_effect.py NCSCLC_comb_input.txt \
            --adj-respB $time $response \
            --N 50000 \
            --predict-type theor \
            --out-prefix ../analysis/NSCLC/Chemo_Atezo_${time}_${response}_theor \
            --extension pdf \
            --out-table \
            --time-max 25 \
            --figsize 6 4
    done
done
