# -*- coding: utf-8 -*-
import h5py
import numpy as np

tmp = h5py.File('zs11.500K_ICE.h5','r') 
np.savetxt('tmp3',tmp['matrix'],fmt='%s', delimiter='\t')
ori_mat=np.loadtxt('tmp3')
#ori_mat=tmp['matrix']
ori_mat_out=[]
for i in range(len(ori_mat)):
        for j in range(len(ori_mat)):
                if (ori_mat[i,j]!=0):
                        ori_mat_out.append(str(i+1)+'\t'+str(j+1)+'\t'+str(ori_mat[i,j]))
                else:
                        pass
np.savetxt('zs11.500K_iced.triple',ori_mat_out,fmt='%s', delimiter='\t')
