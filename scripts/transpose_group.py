#!/usr/bin/env python

"""Transpose group.csv and output sample and group to new file"""

import csv
import sys

import numpy as np
import pandas as pd

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
        out_file = sys.argv[2]
    except:
        print('Please input a file!')
        sys.exit()

    df = pd.read_csv(filename, header=None)
    t_df = np.transpose(df)
    t_df = t_df.drop(columns=[1])
    t_df.columns = ['sample-id', 'group']
    t_df['group'] = t_df['group'].apply(lambda x: 'C' + str(x))
    
    t_df.to_csv(out_file, header=True, sep='\t', index=False)
