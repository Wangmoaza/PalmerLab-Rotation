#!/usr/bin/env bash

source /nas/longleaf/apps/anaconda/2019.10/etc/profile.d/conda.sh
source activate clonetracer

HOME_DIR="/nas/longleaf/home/hhwangbo/Palmer_rotation"
TOOL_DIR="$HOME_DIR/tools/clonetracer_analyze_12"
OUT_DIR="$HOME_DIR/analysis/clonetracer/sample_level"
IN_DIR="$HOME_DIR/data/clonetracer"

files=$(echo `find $IN_DIR/$1 -type f` | tr " " ",")

mkdir -p $OUT_DIR/$1
python $TOOL_DIR/clonTracer_countBarcodes.py -i $files -o $OUT_DIR/${1}/$1
