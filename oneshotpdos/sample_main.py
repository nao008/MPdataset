import  makepdos as mp
import re,os,shutil
from readinfo import setcifdata
import pdosplot as pdosp

#parsing ciffiles in the direcoty
#there result direcoty
print(os.getcwd)
cifdir='/samplecif'
#cifdir='~/sampletcif'

#fs.make_pdos_in_di(cifdir)
mp.listup_cif(cifdir)
ciflist=[s.replace('\n','') for s in open(cifdir+'/cif_list.txt').readlines()]


#saving pdos
try:
    os.mkdir(cifdir+'/pdospng')
except:
    shutil.rmtree(cifdir+'/pdospng')
    os.mkdir(cifdir+'/pdospng')
os.chdir(cifdir+'/pdospng')
cifdatas=[setcifdata(s) for s in ciflist]

for cifdata in cifdatas:
    if cifdata.allsite!=None:
        os.mkdir(cifdata.cifnumber)
        os.chdir(cifdata.cifnumber)
        for j in range(len(cifdata.allsite)):
            figname=str(cifdata.formular)+str(cifdata.allsite[j])
            key='dos.isp1.site{0:03d}.tmp'.format(j+1)
            for i in range(1,4):
                pdosp.savepdos(cifdir+'/result'+'/'+cifdata.cifnumber,orbital=[i],fig_name=figname,key=key)
        os.chdir('..')


import re
import os,shutil
from readinfo import setcifdata
from pdosfilter import pdosfilter as pf
#parsing ciffiles in the direcoty
#there result direcoty
cifdir='/home/fujikazuki/gaustest'
ciflist=[s.replace('\n','') for s in open(cifdir+'/cif_list.txt').readlines()]
cifdatas=[setcifdata(s) for s in ciflist]

try:
    os.mkdir(cifdir+'/p_orbital_pdos')
except:
    shutil.rmtree(cifdir+'/p_orbital_pdos')
    os.mkdir(cifdir+'/p_orbital_pdos')
os.chdir(cifdir+'/p_orbital_pdos')

"""
cifdatas=[setcifdata(s) for s in ciflist]
for cifdata in cifdatas:
    if cifdata.allsite!=None:
        os.mkdir(cifdata.cifnumber)
        os.chdir(cifdata.cifnumber)
        pdosdata=pdosfilter(cifdata.resultadress)
        for j,s in enumerate(cifdata.allsite):
            figname=cifdata.formular+str(s)
            key='dos.isp1.site{0:03d}.tmp'.format(j+1)
            for i in range(1,4):
                pdosdata.savepdos(orbital=[i],fig_name=figname,key=key)
        os.chdir('..')
"""
import numpy as np

for sigma in np.arange(0.05,1,0.05):
    
    print('sigma=',sigma)
    try:
        os.mkdir("sigma="+str(sigma))
    except:
        shutil.rmtree("sigma="+str(sigma))
        os.mkdir("sigma="+str(sigma))
    os.chdir("sigma="+str(sigma))
    for cifdata in cifdatas:
        print(cifdata.cifnumber)
        if cifdata.allsite!=None:
            try:
                os.mkdir(cifdata.cifnumber)
            except:
                shutil.rmtree(cifdata.cifnumber)
                os.mkdir(cifdata.cifnumber)
            os.chdir(cifdata.cifnumber)
            pdosdata=pf(cifdata.resultadress)
            pdosdata.gaussianfilter(sigma=sigma)
            title=cifdata.formular
            pdosdata.savepdos_sameorbital(cifdata.specsite,fig_name=title,orbital=[1,2,3],outputcsv=True)
            os.chdir('..')
    os.chdir('..')
    break