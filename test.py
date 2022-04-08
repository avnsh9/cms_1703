import pandas as pd
data=pd.read_csv('thetas.dat',delim_whitespace=True,header=None)
d=data.sort_values(by=0)
print(d)
d.to_csv('thetass.dat',sep=' ',index=False,header=False)