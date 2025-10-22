# -*- coding: utf-8 -*-
import os
import sys
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Name')
parser.add_argument('--name', '-n', help='Name', required=True)
args = parser.parse_args()
name=args.name

tmp = pd.read_csv("1."+name,sep="\t",header=None)
tmp = np.array(tmp).tolist()

finalBlock = []
FlagChr = tmp[0][0]
FlagStart = tmp[0][1]
FlagEnd = tmp[0][2]

for i in range(1,len(tmp)) :
    newchr = tmp[i][0]
    newstart = tmp[i][1]
    newend = tmp[i][2]
    if ((newchr == FlagChr) & (newstart > FlagEnd)) | (newchr != FlagChr) : 
        t0 = FlagChr
        t1 = FlagStart
        t2 = FlagEnd
        finalBlock.append([t0,t1,t2])
        FlagStart = newstart
        FlagEnd = newend
    elif (newchr == FlagChr) & (newend > FlagEnd):
        FlagEnd =  newend
    if i == len(tmp)-1 :
        t0 = FlagChr
        t1 = FlagStart
        t2 = FlagEnd
        finalBlock.append([t0,t1,t2])

fr = open("2."+name,"w")
for i in range(len(finalBlock)):
    fr.write( finalBlock[i][0]+'\t'+str(finalBlock[i][1])+'\t'+str(finalBlock[i][2])+'\n' )

