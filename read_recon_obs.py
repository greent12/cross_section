import pandas as pd

def read_file(file):
   names = ['var','ID','type','time','lat','lon','pressure','usag','val1','inc1','val2','inc2']
   data = pd.read_csv(file,usecols=[0,2,4,5,6,7,8,9,10,11,12,13],header=None,names=names,delim_whitespace=True)
   return data.time.values.tolist(),data.lat.values.tolist(),data.lon.values.tolist(),data.pressure.values.tolist()


