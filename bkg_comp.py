from turtle import color
import numpy as np
import matplotlib.pyplot as plt

calculated_monoj1=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=17,skip_footer=7)
calculated_monoj2=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=14,skip_footer=26)
calculated_monov=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=36)
calculated_monoj=np.concatenate((calculated_monoj1,calculated_monoj2))
a=np.sum(calculated_monoj)
b=np.sum(calculated_monov)
print("sum of calculated events in all bins")
print(a+b)
calculated_monoj1_err=np.genfromtxt("/home/avnsh9/workspace/checkmate2/results/My_New_Run/evaluation/myprocess_processResults.txt",usecols=3,skip_header=4,skip_footer=7)
calculated_monoj2_err=np.genfromtxt("/home/avnsh9/workspace/checkmate2/results/My_New_Run/evaluation/myprocess_processResults.txt",usecols=3,skip_header=1,skip_footer=26)
calculated_monov_err=np.genfromtxt("/home/avnsh9/workspace/checkmate2/results/My_New_Run/evaluation/myprocess_processResults.txt",usecols=3,skip_header=23)
calculated_monoj_err=np.concatenate((calculated_monoj1_err,calculated_monoj2_err))
print(calculated_monoj_err)
print(calculated_monov_err)



given_monoj=np.genfromtxt('z_muon_given.dat',usecols=2,skip_header=1)
given_monov=np.genfromtxt('z_muon_given.dat',usecols=6,skip_footer=15,skip_header=1)
c=np.sum(given_monoj)
d=np.sum(given_monov)
print("sum of given events in all bins")
print(c+d)


x=[i+1 for i in range(22)]


x_monoj=np.genfromtxt('z_muon_given.dat',usecols=3,skip_header=1)


plt.bar(x,given_monoj)
plt.bar(x,calculated_monoj,width=0.4,yerr=calculated_monoj_err)

plt.legend(['given','calculated'])
plt.xlabel('Hadronic_recoil bins no: 1 to 22')
plt.ylabel('#events')
plt.title('comparison of events in monojet bins')





plt.show()

x_v=[i+1 for i in range(7)]
plt.bar(x_v,given_monov)
plt.bar(x_v,calculated_monov,width=0.4,yerr=calculated_monov_err)
plt.legend(['given','calculated'])
plt.xlabel('hadronic_recoil bins : 1 to 7')
plt.ylabel('#events')
plt.title('comparison of events in mono_V bins')
plt.show()

ratio_j=calculated_monoj/given_monoj
ratio_v=calculated_monov/given_monov

ratio_j=np.round(ratio_j,decimals=2)
ratio_v=np.round(ratio_v,decimals=2)

plt.plot(x,ratio_j)
plt.axhline(y=1, color='r', linestyle='-')

plt.xlabel('hadronic_recoil bins : 1 to 22')
plt.ylabel('calculated/given')
plt.title('ratio in monojet bins')
for a,b in zip(x,ratio_j):
    plt.text(a,b,b)
plt.show()


plt.plot(x_v,ratio_v)
plt.axhline(y=1, color='r', linestyle='-')
plt.xlabel('hadronic_recoil bins : 1 to 7')
plt.title('ratio in mono_V bins')
plt.ylabel('calculated/given')
for a,b in zip(x_v,ratio_v):
    plt.text(a,b,b)
plt.show()

# scatter plot

had_recoil=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/analysisstdout_cms1703_test.log',usecols=1)
muons_pt=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/analysisstdout_cms1703_test.log',usecols=0)
plt.scatter(had_recoil,muons_pt)

plt.ylabel('muons_pt')
plt.xlabel('hadronic_recoil')
plt.title('scatter plot')
plt.axis([100,1000,100,1000])
plt.show()



