import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
lakh_eleven=56.89
calculated_monoj1=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=27,skip_footer=7)
calculated_monoj2=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=14,skip_footer=26)
calculated_monoj=np.concatenate((calculated_monoj1,calculated_monoj2))
norm=lakh_eleven/calculated_monoj[0]

calculated_monoj=calculated_monoj*norm
given_monoj=np.genfromtxt('z_muon_given.dat',usecols=2,skip_header=11)
ratio=calculated_monoj/given_monoj
print(ratio)

