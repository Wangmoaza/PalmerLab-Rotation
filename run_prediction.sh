#!/usr/bin/env bash

python predict_comb_therapy_effect.py \
    --N 50000 \
    --predict-type theor \
    --out-prefix ../analysis/NSCLC/Chemo_Atezo_original_theor \
    --extension pdf \
    --out-table \
    --time-max 25 \
    --figsize 6 4 \
    NCSCLC_comb_input.txt
