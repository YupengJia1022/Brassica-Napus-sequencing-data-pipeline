# -*- coding: utf-8 -*-

import matplotlib
import numpy as np
import h5py,os,argparse
import pandas as pd
import matplotlib.pyplot as plt
import iced

# 补集
def buji(l,i):
    a=[]
    for j in l:
        if j != i:
            a.append(j)
    return(a)

# 交集
def jiaoji(l,m):
    a=[]
    for i in l:
        if i in m:
            a.append(i)
    return(a)

# CPMG标准化
def cpmg(matrix, resolution):
    genomeSize = matrix.shape[0] * resolution / 1e9
    allValidPairs = np.nansum(matrix) / 1e6
    c = allValidPairs / genomeSize
    newMatrix = matrix / c
    return(newMatrix)

Args=argparse.ArgumentParser(description='juicer out matrix to a complete matrix')
Args.add_argument('-i','--indir',dest='indir',help='juicer output matrix path')
Args.add_argument('-b','--bed',dest='bed',help='a bed-format file of bin information')
Args.add_argument('-r','--reso',dest='reso',help='resolution',type=int)
Args.add_argument('-o','--out',dest='out',help='output file prefix')
args=Args.parse_args()

bed=args.bed
indir=args.indir
reso=args.reso
outfile=args.out
bedinfo=pd.read_table(bed,names=['chr','start','end','binid'])
allchrs=list(bedinfo['chr'])
chrs=sorted(list(bedinfo['chr'].drop_duplicates()))
start_i=0
starts={}
for i in chrs:
    starts.setdefault(i,start_i)
    start_i+=allchrs.count(i)
starts.setdefault(i,start_i)
k=0

for chr_i in chrs:
    contact_file = os.path.join(indir,"{0}_{0}.txt".format(chr_i))
    dat=np.loadtxt(contact_file)
    dat[:,0]=np.int64(dat[:,0]/reso)+starts[chr_i]
    dat[:,1]=np.int64(dat[:,1]/reso)+starts[chr_i]
    chrs_others=buji(chrs,chr_i)
    ii=chrs.index(chr_i)
    for chr_j in chrs_others:
        jj=chrs.index(chr_j)
        if jj>=ii:
            contact_file = os.path.join(indir,"{0}_{1}.txt".format(chr_i,chr_j))
            dat2=np.loadtxt(contact_file)
            dat2[:,0]=np.int64(dat2[:,0]/reso)+starts[chr_i]
            dat2[:,1]=np.int64(dat2[:,1]/reso)+starts[chr_j]
        else:
            contact_file = os.path.join(indir,"{1}_{0}.txt".format(chr_i,chr_j))
            dat2=np.loadtxt(contact_file)
            dat2[:,0]=np.int64(dat2[:,0]/reso)+starts[chr_j]
            dat2[:,1]=np.int64(dat2[:,1]/reso)+starts[chr_i]
        dat=np.vstack([dat,dat2])
    if k==0:
        res=dat
        k+=1
    else:
        res=np.vstack([res,dat])
res2=np.array([res[:,1],res[:,0],res[:,2]])
res2=res2.transpose()
res=np.vstack([res,res2])

size=int(res[:,:2].max())+1
matrix=np.zeros((size,size))
x=np.int64(res[:,0])
x_max=np.max(x)
x2=x_max-x
y=np.int64(res[:,1])
matrix[x,y]=res[:,2]
start_sites=[starts[i] for i in chrs]

# 原始交互数据，hdf5文件，包含染色体号、每个染色体在大矩阵里面的起始位置、全基因组的交互矩阵
out1=h5py.File('{}_NONE.h5'.format(outfile),'w')
out1.create_dataset(name='start',data=start_sites)
out1.create_dataset(name='chr',data=chrs)
out1.create_dataset(name='matrix',data=matrix)
out1.close()

# CPMG标准化的交互数据，hdf5文件，包含染色体号、每个染色体在大矩阵里面的起始位置、全基因组的CPMG标准化交互矩阵
cpmgMatrix = cpmg(matrix, reso)
out2=h5py.File('{}_CPMG.h5'.format(outfile),'w')
out2.create_dataset(name='chr', data=chrs)
out2.create_dataset(name='start', data=start_sites)
out2.create_dataset(name='matrix', data=cpmgMatrix)
out2.close()

# ICE标准化的交互数据，hdf5文件，包含染色体号、每个染色体在大矩阵里面的起始位置、全基因组的ICE标准化交互矩阵
iceMatrix = iced.normalization.ICE_normalization(matrix)
out3=h5py.File('{}_ICE.h5'.format(outfile),'w')
out3.create_dataset(name='chr', data=chrs)
out3.create_dataset(name='start', data=start_sites)
out3.create_dataset(name='matrix', data=iceMatrix)
out3.close()
