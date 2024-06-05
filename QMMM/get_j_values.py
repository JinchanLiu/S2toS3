#!/usr/bin/env python3
import numpy as np
np.set_printoptions(suppress=True)

# number of unpaired electrons on each fragment
sys_spin = [3,3,4,3]
# electron configuration, a=spin up ; b=spin down
sys_conf = ['aaaa','aaab','aaba','abaa','baaa','abab','baab','bbaa']
# E(BS)-E(HS), kcal/mol, relative energy, HS not included
rel_e = [-181.4209757,-576.5869190,-135.7846891,164.8264565,-331.8471447,-52.5102651,-345.5639097]

# build the spin configuration lists
def convert_spin(spin,conf):
    sys=[]
    for c in conf:
        sys.append([])
        if len(spin)!=len(c):
            print(spin,c)
            raise ValueError("spinnum and sysconf not match")
        for i in range(len(spin)):
            sys[-1].append(spin[i])
            if c[i]=='b':
                sys[-1][i]*=-1
    print(sys)
    return sys

def get_J_values(unpair_ele,relative_energy):
    n_center=len(unpair_ele[0])
    n_bs=2**(n_center-1)-1
    n_pairs=n_center*(n_center-1)//2
    print(n_center,n_bs,n_pairs)
    
    pairs=[]
    for i in range(len(unpair_ele)):
        tmp=[]
        for j in range(n_center-1):
            for k in range(j+1,n_center):
                tmp.append(unpair_ele[i][j]*unpair_ele[i][k]/4.)
        pairs.append(tmp)
    print(pairs)
    mat=np.zeros((n_bs,n_pairs))
    for i in range(1,len(unpair_ele)):
        mat[i-1]=np.array(pairs[i])-np.array(pairs[0])
    #print(mat)
    mat*=-2
    
    u,s,vh=np.linalg.svd(mat)
    #print(u,'\n',s,'\n',vh)
    #print(u.shape,s.shape,vh.shape)
    #print(s,1/s)
    s_mat=np.zeros((n_bs,n_pairs))
    np.fill_diagonal(s_mat,s)
    #print(s_mat)
    #print(u@s_mat@vh)
    
    s_mat_rev=np.zeros((n_pairs,n_bs))
    np.fill_diagonal(s_mat_rev,1./s)
    result1=vh.T@s_mat_rev@u.T@np.array(relative_energy)
    print(result1)      # J values
    print(mat@result1)  # should reproduce the rel_e
    
    # a different way of SVD, for verification
    u2,s2,vh2=np.linalg.svd(mat,full_matrices=False)
    #print(u2,'\n',s2,'\n',vh2)
    #print(u2.shape,s2.shape,vh2.shape)
    #print(u2@np.diag(s2)@vh2)
    result2=vh2.T@np.diag(1./s2)@u2.T@np.array(relative_energy)
    print(result2)      # J values
    print(mat@result2)  # should reproduce the rel_e
    
    #return result1

sys=convert_spin(sys_spin,sys_conf)
get_J_values(sys,rel_e)

# J constants with all 16 configurations
sys_spin2=[3,3,4,3,1]
sys_conf2=['aaaaa','aaaba','aabaa','abaaa','ababa','baaaa','baaba','bbaaa','aaaab','aaabb','aabab','abaab','ababb','baaab','baabb','bbaab']
rel_e2=[-181.284462,-576.572028,-135.815981,-331.907824, 164.727146, -52.199636,-345.690604,  -0.363612,-181.784588,-576.950531,-136.148301,-332.210757, 164.462844, -52.873877,-345.927522]
sys2=convert_spin(sys_spin2,sys_conf2)
get_J_values(sys2,rel_e2)

# J constants with all 16 configurations, very tight
sys_spin2=[3,3,4,3,1]
sys_conf2=['aaaaa','aaaba','aabaa','abaaa','ababa','baaaa','baaba','bbaaa','aaaab','aaabb','aabab','abaab','ababb','baaab','baabb','bbaab']
rel_e2=[-181.652463,-576.584619,-135.782995,-332.254332, 164.791355, -52.568423,-345.550925,   0.003139,-181.424787,-576.596809,-135.816805,-332.029409, 164.811394, -52.585050,-345.551266]
sys2=convert_spin(sys_spin2,sys_conf2)
get_J_values(sys2,rel_e2)

