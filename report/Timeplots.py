"""
get from output files
"""

with open('../src/results/sing-output.dat','r') as f:
	data=f.read()	

data=data.split('\n')

catdat=[]
for i in range(192/2):
	catdat.append(data[2*i]+' '+data[2*i+1])

actdata=[]
for j in catdat:
	actdata.append([int(j.split(' ')[4]),float(j.split(' ')[7][:-1])])	

cnt=-1
sorteddata={0:[],1:[],2:[],3:[]}
for i in actdata:
	if i[0]==1:
		cnt+=1
	sorteddata[cnt].append(i)

import numpy as np
from scipy.optimize import curve_fit
actdata.sort()
xfull=[]
yfull=[]
for i in actdata:
	xfull.append(i[0])
	yfull.append(i[1])

#def func(x,a):
#	return a * np.exp(-b*x) + c
#	#return a / (x)

#popt, pcov = curve_fit(func, xfull, yfull)

def returnvals(actdat):
	x=[]
	y=[]
	for i in actdat:
		x.append(i[0])
		y.append(i[1])
	return x,y

import matplotlib.pyplot as plt
plt.figure(figsize=(7,5))
for j in range(4):
	xt,yt=returnvals(sorteddata[j])
	plt.plot(xt,yt,'o',label='Job '+str(j+1))

#plt.plot(xfull, func(np.array(xfull), *popt), 'k-', label="Fitted Curve")

plt.legend(loc=9)
plt.ylabel('Time (seconds)',size=20)
plt.xlabel('Threads',size=20)
plt.title('Performance of single node (small data set)')
plt.savefig('SingleTime.pdf')


##now do multi file
with open('../src/results/mult-output.dat','r') as f:
	data=f.read()	

data=data.split('\n')

catdat=[]
for i in range(10/2):
	catdat.append(data[2*i]+' '+data[2*i+1])

actdata=[]
for j in catdat:
	actdata.append([int(j.split(' ')[4]),float(j.split(' ')[7][:-1])])

xfull=[]
yfull=[]
for i in actdata:
	xfull.append(i[0])
	yfull.append(i[1])

plt.clf()
plt.figure(figsize=(7,5))
x = np.r_[0:1:25, 90:10:200]
fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)


for j in range(4):
	xt,yt=returnvals(sorteddata[j])
	ax.plot(xt,yt,'o',label='Job '+str(j+1)+' Single Node')

ax2.plot(xfull,yfull,'ko',label='8 Node Perf')
ax2.plot(xfull,yfull,'k')

ax.set_xlim(0,25)
ax.legend(loc=9)
ax2.set_xlim(90,200)
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.yaxis.tick_right()
plt.subplots_adjust(wspace=0.15)

plt.legend(loc=9)
#plt.ylabel('Time (seconds)',size=20)
#plt.xlabel('Threads',size=20)
plt.suptitle('Combined performance (small data set)',size=24)
fig.text(0.5, 0.04, 'Threads', ha='center',size=20)
fig.text(0.04, 0.5, 'Time (seconds)', va='center', rotation='vertical',size=20)
plt.savefig('CombineTimes.pdf')


plt.clf()
plt.plot(xfull,yfull,'ko')
plt.plot(xfull,yfull,'k')
plt.ylabel('Time (seconds)',size=20)
plt.xlabel('Threads',size=20)
plt.title('Performance over 8 nodes (small data set)')
plt.savefig('MultTimes.pdf')

