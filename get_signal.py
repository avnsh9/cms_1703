import numpy as np
calculated_monoj1=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/dark/analysis/myprocess_cms1703_signal.dat',usecols=4,skip_header=11,skip_footer=28)
calculated_monoj2=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/dark/analysis/myprocess_cms1703_signal.dat',usecols=4,skip_header=22,skip_footer=17)
calculated_monoj3=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/dark/analysis/myprocess_cms1703_signal.dat',usecols=4,skip_header=26,skip_footer=7)
calculated_monoj4=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/dark/analysis/myprocess_cms1703_signal.dat',usecols=4,skip_header=12,skip_footer=18)
calculated_monoj5=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/dark/analysis/myprocess_cms1703_signal.dat',usecols=4,skip_header=23,skip_footer=14)
print(calculated_monoj1)
print(calculated_monoj2)
print(calculated_monoj3)
print(calculated_monoj4)
print(calculated_monoj5)

# calculated_monoj=np.concatenate((calculated_monoj1,calculated_monoj2,calculated_monoj3,calculated_monoj4,calculated_monoj5))
# print(calculated_monoj)
calculated_monoj=np.concatenate((calculated_monoj3,calculated_monoj4,calculated_monoj5))
print(calculated_monoj)
a=np.array([152.898, 110.269, 82.4174,  83.5542 , 61.3867,  64.2287,  48.3136,  31.2618 , 25.5778 , 16.4835
, 21.0306 , 14.7783 , 11.3679  ,10.2311 ,  2.84198,  2.27358 , 1.70519 , 3.97877,
  0.    ,   2.84198,  1.70519 , 4.54717])
# print(a)
# print(a[21])
