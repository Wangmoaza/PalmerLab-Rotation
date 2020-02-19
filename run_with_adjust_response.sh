#!/usr/bin/env bash

for time in 1.5 3 4.5
do
    for response in 60 70 80 90
    do
        python predict_comb_therapy_effect.py NCSCLC_comb_input.txt \
            --adj-respB $time $response \
            --N 50000 --predict-type theor \
            --out-prefix ../analysis/NSCLC/Chemo_Atezo_${time}_${response}_theor \
            --extension png
    done
done
