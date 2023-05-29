import glob
import pandas as pd
import tqdm
import re

datalist = glob.glob('MP-DB/mp-*')

dict = {}
for data in tqdm.tqdm(datalist):
    # print(data)
    df = pd.read_csv(data,index_col=0)

    index = df.index

    index_first = index[0]
    index_last = index[-1]

    index_range = (index_first, index_last)
    
    dataname = data.replace('MP-DB/','').replace('.csv','')
    if index_range in dict.keys():
        dict[index_range].append(dataname)
    else:
        dict[index_range] = [dataname]
    


print(dict.keys())