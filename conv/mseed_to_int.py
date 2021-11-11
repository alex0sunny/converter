import os
from obspy import *
import numpy as np

# original mseed files dir
dir = 'd:/converter_data/temp'	# !! use /, not \ !!
k = 1000.0

files = [x for x in os.listdir(dir) if x.endswith('.mseed')]
os.chdir(dir)
for file in files:
    print(file)
    st = read(file)
    for tr in st:
      tr.data = (tr.data * k).astype('int32')
    st.write(file)