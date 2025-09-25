#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
from pathlib import Path

def calculate_methylation_stats(infile, outfile):
    """计算每个bin的甲基化统计"""
    
    # 读取数据
    df = pd.read_csv(infile, sep='\t', header=None,
                    names=['Chr', 'Start', 'End', 'bin_id',
                           'meth_chr', 'meth_pos', 'meth_pos_end', 'strand',
                           'methy_C', 'non_methy_C', 'all_C', 'methy_rate',
                           'Context', 'others'])
    
    # 转换数据类型
    df['methy_C'] = pd.to_numeric(df['methy_C'], errors='coerce').fillna(0)
    df['all_C'] = pd.to_numeric(df['all_C'], errors='coerce').fillna(0)
    
    # 按bin分组计算统计
    results = []
    for bin_id, group in df.groupby('bin_id'):
        bin_info = group.iloc[0]
        
        # 初始化统计
        stats = {
            'CG': {'sites': 0, 'methy_ratio': 'nan'},
            'CHG': {'sites': 0, 'methy_ratio': 'nan'},
            'CHH': {'sites': 0, 'methy_ratio': 'nan'}
        }
        
        # 计算各context类型的统计
        for context in stats.keys():
            context_data = group[group['Context'] == context]
            site_count = len(context_data)
            
            if site_count > 0:
                total_methy = context_data['methy_C'].sum()
                total_reads = context_data['all_C'].sum()
                
                if total_reads > 0:
                    methy_ratio = total_methy / total_reads
                    stats[context]['methy_ratio'] = f"{methy_ratio:.4f}"
                else:
                    stats[context]['methy_ratio'] = 'nan'
            else:
                stats[context]['methy_ratio'] = 'nan'
            
            stats[context]['sites'] = site_count
        
        # 构建结果行
        result_line = [
            bin_info['Chr'], bin_info['Start'], bin_info['End'],
            stats['CG']['sites'], stats['CHG']['sites'], stats['CHH']['sites'],
            stats['CG']['methy_ratio'], stats['CHG']['methy_ratio'], stats['CHH']['methy_ratio']
        ]
        results.append(result_line)
    
    # 保存结果
    with open(outfile, 'w') as f:
        for result in results:
            f.write('\t'.join(map(str, result)) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Calculate methylation statistics')
    parser.add_argument('--name', '-n', required=True, help='Chromosome name')
    args = parser.parse_args()
    
    infile = f"3.{args.name}"
    outfile = f"4.{args.name}"
    
    if Path(infile).exists():
        calculate_methylation_stats(infile, outfile)
        print(f"Processed {args.name}: {outfile}")
    else:
        print(f"Error: Input file {infile} not found")

if __name__ == '__main__':
    main()
