from oneshotpdos.readinfo import set_pdosdata,read_pdos
import matplotlib.pyplot as plt
import os
import glob
import re
import pandas as pd
import tqdm

mpdatalist = glob.glob('../mp-*')

s_count = 0
for mpdata in tqdm.tqdm(mpdatalist):
    # print(mpdata)
    mpnum = re.findall(r'[0-9]+', mpdata)[0]

    sitedatalist = glob.glob(f'{mpdata}/dos.isp1.site*')

    for sitedata in sitedatalist:
        if '.eps' in sitedata:
            print(f'eps file: {sitedata}')
            continue
        s_count += 1
        sitenum = re.search(r'site(\d+)', sitedata).group(1)
        pdosdf=read_pdos(f'{sitedata}')

        pdosdf.to_csv(f'MP-DB/mp-{mpnum}_site{sitenum}.csv')

print(f'pdos data count: {s_count}')


