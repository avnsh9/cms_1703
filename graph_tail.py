from turtle import color
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

lakh_eleven=56.89

calculated_monoj1=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/lakh_pit/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=17,skip_footer=16)
calculated_monoj11=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=27,skip_footer=7)
calculated_monoj2=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=14,skip_footer=26)
calculated_monoj=np.concatenate((calculated_monoj11,calculated_monoj2))
norm=lakh_eleven/calculated_monoj[0]
calculated_monoj=calculated_monoj*norm
calculated_monoj=np.concatenate((calculated_monoj1,calculated_monoj))
print(calculated_monoj)


calculated_monov=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/lakh_pit/analysis/myprocess_cms1703_test_signal.dat',usecols=4,skip_header=36)

a=np.sum(calculated_monoj)
b=np.sum(calculated_monov)
# print("sum of calculated events in all bins")
# print("monoj",a)
# print("monov",b)
calculated_monoj1_err=np.genfromtxt("/home/avnsh9/workspace/checkmate2/results/lakh_pit/evaluation/myprocess_processResults.txt",usecols=3,skip_header=4,skip_footer=16)
calculated_monoj11_err=np.genfromtxt("/home/avnsh9/workspace/checkmate2/results/My_New_Run/evaluation/myprocess_processResults.txt",usecols=3,skip_header=14,skip_footer=7)
calculated_monoj2_err=np.genfromtxt("/home/avnsh9/workspace/checkmate2/results/My_New_Run/evaluation/myprocess_processResults.txt",usecols=3,skip_header=1,skip_footer=26)
calculated_monov_err=np.genfromtxt("/home/avnsh9/workspace/checkmate2/results/lakh_pit/evaluation/myprocess_processResults.txt",usecols=3,skip_header=23)
calculated_monoj_err=np.concatenate((calculated_monoj11_err,calculated_monoj2_err))
calculated_monoj_err=calculated_monoj_err*norm
calculated_monoj_err=np.concatenate((calculated_monoj1_err,calculated_monoj_err))
print(calculated_monoj_err)




given_monoj=np.genfromtxt('z_muon_given.dat',usecols=2,skip_header=1)
given_monov=np.genfromtxt('z_muon_given.dat',usecols=6,skip_footer=15,skip_header=1)
c=np.sum(given_monoj)
d=np.sum(given_monov)
# print("sum of given events in all bins")
# print("givenj",c)
# print("givenv",d)
g_j_e=given_monoj*0.05
g_v_e=given_monov*0.05


x=[i+1 for i in range(22)]


x_monoj=np.genfromtxt('z_muon_given.dat',usecols=3,skip_header=1)
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2.50

ax1 = plt.subplot()
plt.bar(x,given_monoj,yerr=g_j_e,ecolor='r')
plt.bar(x,calculated_monoj,width=0.4,yerr=calculated_monoj_err)
# ticks = [str(i) for i  in [0, 0, 335, 530, 765, 1055, 1200]]
ticks = [str(i) for i  in [0, 200, 400, 600, 800, 1000, 1200]]
ax1.set_xticklabels(ticks) 
print(calculated_monoj)




plt.legend(['Published','Calculated'])
plt.xlabel('Hadronic_recoil (GeV)')
plt.ylabel('#events')





plt.show()

x_v=[i+1 for i in range(7)]
plt.bar(x_v,given_monov,yerr=g_v_e,ecolor='r')
plt.bar(x_v,calculated_monov,width=0.4,yerr=calculated_monov_err)
plt.legend(['Published','Calculated'])
plt.xlabel('hadronic_recoil bins : 1 to 7')
plt.ylabel('#events')

plt.show()

ratio_j=calculated_monoj/given_monoj
ratio_v=calculated_monov/given_monov

ratio_j=np.round(ratio_j,decimals=2)
ratio_v=np.round(ratio_v,decimals=2)
j_r_e=((calculated_monoj_err/given_monoj)**2 + (0.05*ratio_j)**2)**(1/2)
v_r_e=((calculated_monov_err/given_monov)**2 + (0.05*ratio_v)**2)**0.5
print(j_r_e)
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2.50

ax1 = plt.subplot()
 
plt.errorbar(x,ratio_j,yerr=j_r_e,fmt='o',ecolor='r',markersize=15)
plt.axhline(y=1, color='k', linestyle='--')
# ticks = [str(i) for i  in [0, 0, 335, 530, 765, 1055, 1200]]
ticks = [str(i) for i  in [0, 200, 400, 600, 800, 1000, 1200]]

ax1.set_xticklabels(ticks)

plt.xlabel('Hadronic_recoil (GeV)')
plt.ylabel('Calculated/Published')

# for a,b in zip(x,ratio_j):
#     plt.text(a,b,b)
plt.show()


# plt.errorbar(x_v,ratio_v,yerr=v_r_e,fmt='o',ecolor='r')
# plt.axhline(y=1, color='r', linestyle='-')
# plt.xlabel('hadronic_recoil bins : 1 to 7')
# plt.title('ratio in mono_V bins')
# plt.ylabel('calculated/given')
# for a,b in zip(x_v,ratio_v):
#     plt.text(a,b,b)
# plt.show()

# scatter plot

# had_recoil=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/analysisstdout_cms1703_test.log',usecols=1)
# muons_pt=np.genfromtxt('/home/avnsh9/workspace/checkmate2/results/My_New_Run/analysis/analysisstdout_cms1703_test.log',usecols=0)
# plt.scatter(had_recoil,muons_pt)

# plt.ylabel('muons_pt')
# plt.xlabel('hadronic_recoil')
# plt.title('scatter plot')
# plt.axis([100,1000,100,1000])
# plt.show()




