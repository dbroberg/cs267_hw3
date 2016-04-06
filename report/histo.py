"""
Plot histogram of starting kmers for short and long list
"""

shorttest=  [[1.0000 , 451417.7000 ,   341],
[451417.7000 , 902834.4000 ,   595],
[902834.4000 , 1354251.1000 ,   587],
[1354251.1000 , 1805667.8000 ,   665],
[1805667.8000 , 2257084.5000 ,   546],
[2257084.5000 , 2708501.2000 ,   465],
[2708501.2000 , 3159917.9000 ,   623],
[3159917.9000 , 3611334.6000 ,   623],
[3611334.6000 , 4062751.3000 ,   633],
[4062751.3000 , 4514168.0000 ,   658]]

listset=[]
numset=[]
tot=0
for i in shorttest:
	listset.append(str(int(i[0]))+' - '+str(int(i[1])))
	numset.append(i[2])
	tot+=i[2]

perc=[]
for i in shorttest:
	perc.append(round(100*float(i[2])/float(tot),2))

import matplotlib.pyplot as plt
plt.figure(figsize=(20,10))
plt.subplots_adjust(bottom=0.32)
#x=range(int(shorttest[-1][1]))
x=range(len(shorttest))
plt.bar(x,numset)
plt.xticks(x,listset,rotation=45,size=20)
plt.yticks(size=30)
plt.ylabel('Number of kmers',size=30)
plt.xlabel('Address range',size=30)
#plt.title('Histogram of kmers distributed (small set)',size=40)
plt.grid()
for i in range(len(shorttest)):
	plt.text(i+.05,shorttest[i][2]-40,str(perc[i])+' %',color='w',size=20,verticalalignment='center')


plt.savefig('shorthisto.pdf')

longtest = [[1.0000 , 2699925.0000, 18187],
[2699925.0000 , 5399849.0000 , 22391],
[5399849.0000 , 8099773.0000 , 20642],
[8099773.0000 , 10799697.0000 , 20927],
[10799697.0000 , 13499621.0000 , 22963],
[13499621.0000 , 16199545.0000 , 20172],
[16199545.0000 , 18899469.0000 , 20746],
[18899469.0000 , 21599393.0000 , 20735],
[21599393.0000 , 24299317.0000 , 23472],
[24299317.0000 , 26999241.0000 , 20315]]

listset=[]
numset=[]
tot=0
for i in longtest:
	listset.append(str(int(i[0]))+' - '+str(int(i[1])))
	numset.append(i[2])
	tot+=i[2]

perc=[]
for i in longtest:
	perc.append(round(100*float(i[2])/float(tot),2))

plt.clf()
plt.figure(figsize=(20,10))
plt.subplots_adjust(bottom=0.32)
#x=range(int(shorttest[-1][1]))
x=range(len(longtest))
plt.bar(x,numset)
plt.xticks(x,listset,rotation=45,size=20)
plt.yticks(size=30)
plt.ylabel('Number of kmers',size=30)
plt.xlabel('Address range',size=30)
#plt.title('Histogram of kmers distributed (large set)',size=40)
plt.grid()
for i in range(len(longtest)):
	plt.text(i+.38,longtest[i][2]-1000,str(perc[i])+' %',color='w',size=20,verticalalignment='center',horizontalalignment='center')

plt.savefig('longhisto.pdf')


