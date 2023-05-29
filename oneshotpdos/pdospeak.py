from pdosfilter import MpPdosdata as mpd
dir='/home/fujikazuki/Documents/mp2534'

t=mpd(dir)
t.gaussianfilter(sigma=0.5)
t.make_histogram()

print(t.histogram)