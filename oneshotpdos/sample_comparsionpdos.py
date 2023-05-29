from comparsion import ComparsionPdos

dir='/home/fujikazuki/gaustest'
resultdir='/home/fujikazuki/gaustest/classtest'
ciflist=[s.replace('\n','')  for s in open(dir+'/cif_list.txt')]
tc=ComparsionPdos(cifadress_1=ciflist[11],cifadress_2=ciflist[2])
tc.gaussian(Sigma=0.5)
tc.peak_range()
tc.make_L()
import pandas as pd
idx=pd.IndexSlice['s',:]
col=pd.IndexSlice['site_0',:]
import pulp as p
m=p.LpProblem()
#set variable dataframe
val=pd.DataFrame()
from itertools import product
for i,j in  product(tc.L_df.loc[idx,col].droplevel(level=0).index,tc.L_df.loc[idx,col].droplevel(level=0,axis=1).columns):
    val.loc[i,j]=p.LpVariable('w_%i,%i'%(i,j),lowBound=0)


m+=p.lpDot(tc.L_df.loc[idx,col].droplevel(level=0).droplevel(level=0,axis=1),val)

for i in val.index:
    m+=p.lpSum(val.loc[i,:])==1.0
for i in val.columns:
    m+=p.lpSum(val.loc[:,i])==1.0

m.solve()
result_df=pd.DataFrame()

print(val.applymap(p.value))