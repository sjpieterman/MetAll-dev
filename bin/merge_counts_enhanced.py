#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import glob

def merge_featurecounts(files):
    dfs = []
    for f in files:
        sample_name = os.path.basename(f).replace('.featureCounts.txt', '')
        df = pd.read_csv(f, sep='\t', comment='#')
        # featureCounts has columns: Geneid, Chr, Start, End, Strand, Length, Count
        df = df[['Geneid', df.columns[-1]]]
        df.columns = ['Geneid', sample_name]
        df['Geneid'] = df['Geneid'].str.split('.').str[0] # Strip version
        dfs.append(df.set_index('Geneid'))
    return pd.concat(dfs, axis=1)

def merge_star(files):
    dfs = []
    for f in files:
        sample_name = os.path.basename(f).replace('.ReadsPerGene.out.tab', '')
        df = pd.read_csv(f, sep='\t', header=None, skiprows=4)
        # STAR output: 0:gene, 1:unstranded, 2:strand1, 3:strand2
        df = df[[0, 1]]
        df.columns = ['Geneid', sample_name]
        df['Geneid'] = df['Geneid'].str.split('.').str[0]
        dfs.append(df.set_index('Geneid'))
    return pd.concat(dfs, axis=1)

def merge_salmon(files):
    dfs = []
    for f in files:
        # Salmon quant.sf: Name, Length, EffectiveLength, TPM, NumReads
        sample_name = os.path.basename(os.path.dirname(f)).replace('_salmon', '')
        df = pd.read_csv(f, sep='\t')
        df = df[['Name', 'NumReads']]
        df.columns = ['Geneid', sample_name]
        df['Geneid'] = df['Geneid'].str.split('.').str[0]
        dfs.append(df.set_index('Geneid'))
    return pd.concat(dfs, axis=1)

def main():
    parser = argparse.ArgumentParser(description='Merge counts from various tools')
    parser.add_argument('--featurecounts', nargs='*', help='featureCounts output files')
    parser.add_argument('--star', nargs='*', help='STAR ReadsPerGene.out.tab files')
    parser.add_argument('--salmon', nargs='*', help='Salmon quant.sf files')
    parser.add_argument('--output', required=True, help='Output filename')
    
    args = parser.parse_args()
    
    all_dfs = {}
    
    if args.featurecounts:
        all_dfs['featurecounts'] = merge_featurecounts(args.featurecounts)
    if args.star:
        all_dfs['star'] = merge_star(args.star)
    if args.salmon:
        all_dfs['salmon'] = merge_salmon(args.salmon)
        
    if not all_dfs:
        print("No input files provided")
        return

    # Save individual tool matrices
    for tool, df in all_dfs.items():
        tool_out_name = f"{tool}_{args.output}"
        df.to_csv(tool_out_name, sep='\t')
        print(f"Saved {tool} matrix to {tool_out_name}")

    if len(all_dfs) > 1:
        # Create a combined matrix and save as the main output
        combined = pd.concat([df.add_suffix(f'_{tool}') for tool, df in all_dfs.items()], axis=1)
        combined.to_csv(args.output, sep='\t')
        print(f"Saved combined matrix to {args.output}")
    else:
        # Only one tool, save it as the main output (redundant but safe)
        tool = list(all_dfs.keys())[0]
        all_dfs[tool].to_csv(args.output, sep='\t')
        print(f"Saved {tool} matrix as main output: {args.output}")

if __name__ == '__main__':
    main()
