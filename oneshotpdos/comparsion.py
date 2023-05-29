import pandas as pd
from pdosfilter import SameobitalPdosdata as sop
import re,os,itertools
from constant import ORBITAL
from makepdos import listup_cif
import math


def cif_conbination(cif_dir):
    cif_list_txt=cif_dir+'/cif_list.txt'
    if not os.path.isfile(cif_list_txt):
        listup_cif(cif_dir)
    with open(cif_list_txt,mode='r') as f:
            ciflist=[re.search(r"([^/]*?)$",s.strip()).group() for s in f.readlines()]
    return list(itertools.combinations(ciflist,2))

class ComparsionPdos():
    """
    dir='/home/fujikazuki/gaustest'
    resultdir='/home/fujikazuki/gaustest/classtest'
    ciflist=[s.replace('\n','')  for s in open(dir+'/cif_list.txt')]
    tc=ComparsionPdos(cifadress_1=ciflist[11],cifadress_2=ciflist[2])
    tc.gaussian(Sigma=0.5)
    tc.peak_range()
    tc.make_L()
    print(tc.L_df)
    """
    def __init__(self,cifadress_1,cifadress_2):
        self.pdos_data_1=sop(cifadress_1)
        self.pdos_data_2=sop(cifadress_2)
        if len(self.pdos_data_1.specsite)!=len(self.pdos_data_2.specsite):
            print('not same number of specsite')
            return
        self.weight=dict()
        self.specrange=len(self.pdos_data_1.specsite)
    
    def gaussian(self,Sigma):
        self.pdos_data_1.gaussianfilter(sigma=Sigma)
        self.pdos_data_2.gaussianfilter(sigma=Sigma)
    
    def peak_range(self,Ftype='g'):
        self.pdos_data_1.peakdata(ftype=Ftype)
        self.pdos_data_2.peakdata(ftype=Ftype)
    
    def make_L(self):
        siten_dic=dict()
        for siten in [0,1]:
            col_d1=pd.IndexSlice[self.pdos_data_1.atomlist[siten],:]
            col_d2=pd.IndexSlice[self.pdos_data_2.atomlist[siten],:]
            orbital_dic=dict()
            for o in ORBITAL:
                idx=pd.IndexSlice[o,:]
                orbital_df=pd.DataFrame()
                for I,i in self.pdos_data_1.peaks.loc[idx,col_d1].iterrows():
                    for J,j in self.pdos_data_2.peaks.loc[idx,col_d2].iterrows():
                        l=abs(i.droplevel(level=0)-j.droplevel(level=0))
                        L=l['arg_y']*math.sqrt(l['arg_x']*2+l['peak_range']*2)
                        orbital_df.loc[I[1],J[1]]=L
                orbital_dic[o]=orbital_df
            siten_dic['site_%i'%siten]=pd.concat(orbital_dic)
        self.L_df=pd.concat(siten_dic,axis=1)
        
