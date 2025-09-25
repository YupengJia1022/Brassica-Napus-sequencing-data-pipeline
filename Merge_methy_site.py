#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
from pathlib import Path

def merge_adjacent_regions(infile, outfile):
    """合并相邻的低甲基化区域"""
    
    # 读取数据
    df = pd.read_csv(infile, sep='\t', header=None,
                    names=['Chr', 'Start', 'End', 'CG_sites', 'CHG_sites', 'CHH_sites',
                           'CG_methy', 'CHG_methy', 'CHH_methy'])
    
    if len(df) == 0:
        # 如果没有数据，创建空文件
        with open(outfile, 'w') as f:
            pass
        return
    
    # 按染色体和位置排序
    df = df.sort_values(['Chr', 'Start'])
    
    merged_regions = []
    current_chr = df.iloc[0]['Chr']
    current_start = df.iloc[0]['Start']
    current_end = df.iloc[0]['End']
    
    for i in range(1, len(df)):
        row = df.iloc[i]
        
        # 如果相邻且在同一染色体，则合并
        if (row['Chr'] == current_chr and row['Start'] <= current_end):
            current_end = max(current_end, row['End'])
        else:
            # 保存当前区域
            merged_regions.append([current_chr, current_start, current_end])
            # 开始新区域
            current_chr = row['Chr']
            current_start = row['Start']
            current_end = row['End']
    
    # 添加最后一个区域
    merged_regions.append([current_chr, current_start, current_end])
    
    # 保存结果
    with open(outfile, 'w') as f:
        for region in merged_regions:
            f.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")

def main():
    parser = argparse.ArgumentParser(description='Merge adjacent methylation regions')
    parser.add_argument('--name', '-n', required=True, help='Chromosome name')
    args = parser.parse_args()
    
    infile = f"5.{args.name}"
    outfile = f"6.{args.name}"
    
    if Path(infile).exists():
        merge_adjacent_regions(infile, outfile)
        print(f"Merged regions for {args.name}: {outfile}")
    else:
        print(f"Error: Input file {infile} not found")

if __name__ == '__main__':
    main()
