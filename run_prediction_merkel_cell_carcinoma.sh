#!/usr/bin/env bash

python predict_comb_therapy_effect.py \
--out-prefix ~/Dropbox/ExtractIPD/rare_cancers/merkel_cell_carcinoma/chemo_pembrolizumab_combination_theor_kmplot.pdf \
--predict-type theor \
~/Dropbox/ExtractIPD/rare_cancers/merkel_cell_carcinoma/prediction_input.csv

python predict_comb_therapy_effect.py \
--out-prefix ~/Dropbox/ExtractIPD/rare_cancers/merkel_cell_carcinoma/chemo_pembrolizumab_stoch_combination_kmplot.pdf \
--predict-type stoch \
~/Dropbox/ExtractIPD/rare_cancers/merkel_cell_carcinoma/prediction_input.csv
