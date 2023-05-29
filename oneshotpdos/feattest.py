from comparsion import ComparsionPdos 


dir='/home/fujikazuki/gaustest'
resultdir='/home/fujikazuki/gaustest/classtest'
ciflist=[s.replace('\n','')  for s in open(dir+'/cif_list.txt')]
tc=ComparsionPdos(cifadress_1=ciflist[11],cifadress_2=ciflist[2])
tc.gaussian(Sigma=0.5)
tc.peak_range(Ftype='h')
tc.make_L()
print(tc.L_df)