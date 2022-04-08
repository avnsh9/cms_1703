from matplotlib import pyplot as plt
import numpy as np

y=np.loadtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/analysisstdout_cms1703_test.log',usecols=0)
x=np.loadtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/analysisstdout_cms1703_test.log',usecols=1)

#plotting histogram
plt.scatter(x,y)
plt.xlabel('(summing up all the constructed particles and removing lep and photons')
plt.ylabel('(using missingET and removing leptons and photons')
plt.show()
