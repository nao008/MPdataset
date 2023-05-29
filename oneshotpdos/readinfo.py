#read and set cif's infomation from orignal cif file and SiteInfo.lmchk
#set data teyp is dict. format is {cifnumber:[formula structural,{site number:[site atom,position]]}
from scipy.signal import argrelmax,argrelmin
from scipy import integrate
import math
import re,os
import glob
import pandas as pd
from .constant import ORBITAL
import numpy as np

def siteinfo(file):
    linens=[re.sub(' {2,}','/',s.replace('\n','')).split('/') for s in open(file).readlines()]
    
    relist=list()
    for s in linens:
        relist.append([s[3],(s[5],s[6],s[7])])
    return relist

def formularinfo(file):
    testfile=[s.replace('\n','') for s in open(file).readlines()]
    for s in testfile:
        if re.match('_chemical_formula_structural',s):
            return re.sub("_chemical_formula_structural|\n||'| ",'',s)
    for s in testfile:
        if re.match('_chemical_formula_sum',s):
            print("Instead, capture chemical sum in "+re.search(r"([^/]*?)$",file.strip()).group().replace('.cif',''))
            return re.sub("_chemical_formula_|\n||'| ",'',s)

    print('no formula data '+re.search(r"([^/]*?)$",file.strip()).group().replace('.cif',''))
    return None

def specinfo(file):
    #[site number,[siteinfo]]
    linens=[re.sub(' {2,}','/',s.replace('\n','')).split('/') for s in open(file).readlines()]
    relist=dict()
    for s in linens:
        relist[s[2]]=[s[1],[s[3],(s[5],s[6],s[7])]]
    return [s for s in relist.values()]


def listupAllsite(dir):
    '''sample 
    dir='/home/fujikazuki/testd'

    with open(cifdir+'/cif_list.txt',mode='r') as f:
        cifnumber=[re.search(r"([^/]*?)$",s.strip()).group().replace('.cif','') for s in f.readlines()]

    for s in cifnumber:
        print(fs.listupAllsite(cifdir+'/result'+'/'+s))
    '''
    linens=open(dir+'/ciftext.txt').readlines()
    matchlist=list()
    for i,line in enumerate(linens):
        if re.match('All sites',line):
            j=i+2#start point at Atom in All site
            break
    atomicline=list()
    while True:
        if re.match('\n',linens[j]):
            break
        atomicline.append(linens[j].replace('\n',''))
        j+=1
    return atomicline

class setcifdata():
    #only set cif file name 
    #ciffile='~/1000017.cif'
    """
    cifdir='/home/fujikazuki/ciflist_test'
    with open(cifdir+'/cif_list.txt',mode='r') as f:
        ciflist=[s.replace('\n','') for s in f]
    cifdata=[setcifdata(s) for s in ciflist]
    for i in range(len(cifdata)):
        print(cifdata[i].allsite)
    """
    def __init__(self,ciffile):
        self.cifadress=ciffile
        self.cifnumber=re.search(r"([^/]*?)$",ciffile.strip()).group().replace('.cif','')
        self.resultadress=ciffile.replace(self.cifnumber+'.cif','')+'result/'+self.cifnumber
        atomsitefile=self.resultadress+'/SiteInfo.lmchk'
        if os.path.isfile(atomsitefile):
            self.allsite=siteinfo(atomsitefile)
            self.specsite=specinfo(atomsitefile)
        else:
            print('No such file is'+atomsitefile)
            self.allsite=None
        if os.path.isfile(ciffile):
            self.formular=formularinfo(ciffile)
        else:
            print('No such file is '+ciffile)
            self.formular=None


class SetMpData:
    def __init__(self,mpdir):
        self.resultadress=mpdir
        atomsitefile=self.resultadress+'/SiteInfo.lmchk'
        if os.path.isfile(atomsitefile):
            self.allsite=siteinfo(atomsitefile)
            self.specsite=specinfo(atomsitefile)
        else:
            print('No such file is'+atomsitefile)
            self.allsite=None

def read_pdos(file:str) -> pd.DataFrame:
        pdos=pd.read_csv(file,header=None,index_col=0,comment='#',delim_whitespace=True)
        r=13.605
        orbital=ORBITAL
        realpdosdic=pd.DataFrame(columns=orbital,index=pdos.index.to_numpy()*r).fillna(0)
        for i,o in enumerate(orbital,1):
            anl=int(0.5*i*((i-1)*2+2))
            anf=int(0.5*(i-1)*((i-2)*2+2)+1)
            for j in range(anf,anl+1):
                realpdosdic.loc[:,o]+=pdos.loc[:,j].values
        realpdosdic=realpdosdic/r
        return realpdosdic

def read_pdos(file:str) -> pd.DataFrame:
        """sample
        read_pdos(file='mp-149')"""
        pdos=pd.read_csv(file,header=None,index_col=0,comment='#',delim_whitespace=True)
        r=13.605
        orbital=ORBITAL
        realpdosdic=pd.DataFrame(columns=orbital,index=pdos.index.to_numpy()*r).fillna(0)
        for i,o in enumerate(orbital,1):
            anl=int(0.5*i*((i-1)*2+2))
            anf=int(0.5*(i-1)*((i-2)*2+2)+1)
            for j in range(anf,anl+1):
                realpdosdic.loc[:,o]+=pdos.loc[:,j].values
        realpdosdic=realpdosdic/r
        return realpdosdic


def set_pdosdata(directory):
        '''sample
        dir='~/ciflist/result/1528444'
        dict_df=set_pdosdata(dir)
        '''
        if not os.path.isdir(directory):
            print(directory+" dose no exist")
            return
        cifid=os.path.basename(directory)
        files=list()
        for file in glob.glob(directory+'/dos.isp1.site*'):
            if file.endswith(cifid):
                files.append(file)

        pdosdic={key:read_pdos(key) for key in files}
        return pd.concat(pdosdic,axis=0)

def set_sameorbital(specdata,pdosdata):
    """sample specdata=[site number,[site element,(position)]]
    specdata=[['2', ['Cd', ('2.387228', '0.000000', '3.374500')]], ['4', ['S', ('2.387228', '0.000000', '5.972865')]]]
    """
    level_1_keys=np.unique(pdosdata.index.get_level_values(0).to_numpy())
    idx=pd.IndexSlice[level_1_keys,:]
    xdata=pdosdata.xs(level_1_keys[0],level=0).index
    ykeys=[str(s[1]) for s in specdata]
    sameorbital_pdos=dict()
    for o in ORBITAL:
        empty_pdos=pd.DataFrame(index=xdata,columns=ykeys)
        for i,s in enumerate(ykeys):
            sitenumber=int(specdata[i][0])
            idx=pd.IndexSlice['dos.isp1.site{0:03d}.tmp'.format(sitenumber),:]
            empty_pdos[s][:]=pdosdata.loc[idx,o]
        sameorbital_pdos[o]=empty_pdos
    return  pd.concat(sameorbital_pdos)


def peakrange(xdata,ydata):
    argmax=argrelmax(ydata)[0]
    argmin=argrelmin(ydata)[0]
    for i in range(len(argmax)-len(argmin)+1):
        argmax=np.append(argmax,0)
    vii=0
    maxi=0
    resultlist=list()
    for i,v in enumerate(argmin):
        if xdata[vii]<=xdata[argmax[maxi]] and xdata[v]>=xdata[argmax[maxi]]:
            I=integrate.simps(ydata[vii:v],xdata[vii:v])
            #sigma=math.pi*pow(ydata[argmax[maxi]]/I,2)
            sigma=pow(2*math.pi,-0.5)*I/ydata[argmax[maxi]]
            resultlist.append((xdata[maxi],ydata[maxi],sigma*2))
            vii=v
            maxi+=1
        if v==argmin[-1]:
            if xdata[v]<=xdata[argmax[maxi]] and xdata[-1]>=xdata[argmax[maxi]]:
                I=integrate.simps(ydata[v:-1],xdata[v:-1])
                #sigma=math.pi*pow(ydata[argmax[maxi]]/I,2)
                sigma=pow(2*math.pi,-0.5)*I/ydata[argmax[maxi]]
                resultlist.append((xdata[maxi],ydata[maxi],sigma*2))
    return resultlist


def peakrange_histogram(xdata,ydata):
    argmax=argrelmax(ydata)[0]
    argmin=argrelmin(ydata)[0]
    for i in range(len(argmax)-len(argmin)+1):
        argmax=np.append(argmax,0)
    vii=0
    maxi=0
    resultlist=list()
    for i,v in enumerate(argmin):
        if xdata[vii]<=xdata[argmax[maxi]] and xdata[v]>=xdata[argmax[maxi]]:
            I=integrate.simps(ydata[vii:v],xdata[vii:v])
            prange=I/ydata[argmax[maxi]]
            resultlist.append((xdata[maxi],ydata[maxi],prange))
            vii=v
            maxi+=1
        if v==argmin[-1]:
            if xdata[v]<=xdata[argmax[maxi]] and xdata[-1]>=xdata[argmax[maxi]]:
                I=integrate.simps(ydata[v:-1],xdata[v:-1])
                prange=I/ydata[argmax[maxi]]
                resultlist.append((xdata[maxi],ydata[maxi],prange))
    return resultlist