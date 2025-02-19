#!/usr/bin/env python3

import argparse
import csv
import json
import sys
import pandas as pd 


def filter_best_bitscore(df, group_col, score_col):
    idxmax = df.groupby(group_col)[score_col].idxmax()
    df = df.loc[idxmax].reset_index(drop=True)
    return df

def main(args):
    blast_df = pd.read_csv(args.input)
    filtered_df = filter_best_bitscore(blast_df, args.group_col, args.score_col)
    filtered_df.to_csv(args.output, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-g', '--group_col', default='query_seq_id')
    parser.add_argument('-s', '--score_col', default='bitscore')
    
    args = parser.parse_args()
    main(args)
