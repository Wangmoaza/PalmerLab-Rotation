#!/usr/bin/env python

# make experiment level file structure
# move sample folders to PDX_line folders
import os

dirpath = '../../data/clonetracer/'

full_names = [i for i in os.listdir(dirpath) if "GW0" in i]
f = open(dirpath + 'sample_info_simple.tsv', 'r')

f.readline()  # skip header line
for line in f.readlines():
    tokens = line.strip().split('\t')
    sample, pdx = tokens[0], tokens[1]
    for full_name in full_names:
        if sample in full_name:
            os.rename(dirpath + full_name, '{0}{1}/{2}'.format(dirpath, pdx, sample))

f.close()
