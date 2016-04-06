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
#plt.title('Performance of single node (small data set)')
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
#plt.suptitle('Combined performance (small data set)',size=24)
fig.text(0.5, 0.04, 'Threads', ha='center',size=20)
fig.text(0.04, 0.5, 'Time (seconds)', va='center', rotation='vertical',size=20)
plt.savefig('CombineTimes.pdf')

#now redo just the multtimes plot
plt.clf()
plt.plot(xfull,yfull,'ko')
plt.plot(xfull,yfull,'k')
plt.ylabel('Time (seconds)',size=20)
plt.xlabel('Threads',size=20)
#plt.title('Performance over 8 nodes (small data set)')
plt.savefig('MultTimes.pdf')


#now do time scaling for single node data
#get 1node avg
singthread=[]
for i in range(4):
	singthread.append(sorteddata[i][0][1])

singval=np.mean(singthread)

def returnvals(actdat):
	x=[]
	y=[]
	for i in actdat:
		x.append(i[0])
		y.append(100*(float(singval)/i[0])/float(i[1]))
	return x,y

plt.clf()
plt.figure(figsize=(7,5))
for j in range(4):
	xt,yt=returnvals(sorteddata[j])
	plt.plot(xt,yt,'o',label='Job '+str(j+1))

plt.ylim(0,100)
plt.legend(loc=9)
plt.ylabel('Percentage of Ideal Scaling (%)',size=18)
plt.xlabel('Threads',size=18)
#plt.title('Scaling on single node (small data set)',size=21)
plt.savefig('SingleScaling.pdf')





#Now Do similar plots but for BIG data set
with open('../src/results/multBIG-output.dat','r') as f:
	data=f.read()	

data=data.split('\n')

catdat=[]
for i in range(50/2):
	catdat.append(data[2*i]+' '+data[2*i+1])

actdata=[]
for j in catdat:
	actdata.append([int(j.split(' ')[4]),float(j.split(' ')[7][:-1])])	

cnt=-1
sorteddata={0:[],1:[],2:[],3:[],4:[]}
for i in actdata:
	if i[0]==96:
		cnt+=1
	sorteddata[cnt].append(i)

def returnvals(actdat):
	x=[]
	y=[]
	for i in actdat:
		x.append(i[0])
		y.append(i[1])
	return x,y

plt.clf()
plt.figure(figsize=(7,5))
for j in range(5):
	xt,yt=returnvals(sorteddata[j])
	plt.plot(xt,yt,'o',label='Job '+str(j+1))

#plt.plot(xfull, func(np.array(xfull), *popt), 'k-', label="Fitted Curve")

plt.legend(loc=1)
plt.ylabel('Time (seconds)',size=20)
plt.xlabel('Threads',size=20)
#plt.title('Performance of single node (small data set)')
plt.savefig('LargeMultTime.pdf')

#nOw do scaling for BIG data set
with open('../src/results/singleBIG.dat','r') as f:
	data1=f.read()	

data1=data1.split('\n')
singthread=[]
for i in range(3):
	singthread.append(float(data1[2*i+1].split()[-1][:-1]))

singval=np.mean(singthread)  #this is ideal time...
def returnvals(actdat):
	x=[]
	y=[]
	for i in actdat:
		x.append(i[0])
		y.append(100*(float(singval)/i[0])/float(i[1]))
	return x,y

plt.clf()
plt.figure(figsize=(7,5))
for j in range(5):
	xt,yt=returnvals(sorteddata[j])
	plt.plot(xt,yt,'o',label='Job '+str(j+1))

#plt.plot(xfull, func(np.array(xfull), *popt), 'k-', label="Fitted Curve")

plt.legend(loc=1)

plt.ylabel('Percentage of Ideal Scaling (%)',size=18)
plt.xlabel('Threads',size=18)
#plt.title('Performance of single node (small data set)')
plt.savefig('LargeMultScaling.pdf')



#now get the scaling plot for big data on single node...
with open('../src/results/bigset_singlenode.dat','r') as f:
	data=f.read()	

data=data.split('\n')

catdat=[]
for i in range(48/2):
	catdat.append(data[2*i]+' '+data[2*i+1])

actdata1=[]
for j in catdat:
	actdata1.append([int(j.split(' ')[4]),float(j.split(' ')[7][:-1])])	


singval=actdata1[0][1]  #this is ideal time...
def returnvals(actdat):
	x=[]
	y=[]
	for i in actdat:
		x.append(i[0])
		y.append(100*(float(singval)/i[0])/float(i[1]))
	return x,y

plt.clf()
plt.figure(figsize=(7,5))
xt,yt=returnvals(actdata1)
plt.plot(xt,yt,'o')

#plt.plot(xfull, func(np.array(xfull), *popt), 'k-', label="Fitted Curve")

plt.legend(loc=1)

plt.ylabel('Percentage of Ideal Scaling (%)',size=18)
plt.xlabel('Threads',size=18)
#plt.title('Performance of single node (small data set)')
plt.savefig('SingleBIGScaling.pdf')

#now do same plot but with large number of threads included as well?

with open('../src/results/multBIG-output.dat','r') as f:
	data=f.read()	

data=data.split('\n')

catdat=[]
for i in range(50/2):
	catdat.append(data[2*i]+' '+data[2*i+1])

actdata=[]
for j in catdat:
	actdata.append([int(j.split(' ')[4]),float(j.split(' ')[7][:-1])])	

cnt=-1
sorteddata={0:[],1:[],2:[],3:[],4:[]}
for i in actdata:
	if i[0]==96:
		cnt+=1
	sorteddata[cnt].append(i)


x=[]
y=[]
for i in range(len(sorteddata[0])):
	x.append(sorteddata[0][i][0])
	yval=np.mean([sorteddata[0][i][1],sorteddata[1][i][1],sorteddata[2][i][1],sorteddata[3][i][1],sorteddata[4][i][1]])
	y.append(100*(singval/x[-1])/yval)

#def func(x,a,b,c):
#	return a /(c*x-b)

#popt, pcov = curve_fit(func, xfull, yfull)
#xfull=[]
#for i in x:
#	xfull.append(i)


plt.clf()
plt.figure(figsize=(7,5))
xt,yt=returnvals(actdata1)
plt.plot(xt,yt,'o')
plt.plot(x,y,'o')

#for i in xt:
#	xfull.append(i)

#xfull=xfull.sort()
#plt.plot(xfull, func(np.array(xfull), *popt), 'k-', label="Fitted Curve")
plt.legend(loc=1)
plt.ylabel('Percentage of Ideal Scaling (%)',size=18)
plt.xlabel('Threads',size=18)
#plt.title('Performance of single node (small data set)')
plt.savefig('MultiANDSingleBIGScaling.pdf')


#Now with time...

def returnvals(actdat):
	x=[]
	y=[]
	for i in actdat:
		x.append(i[0])
		y.append(float(i[1]))
	return x,y

plt.clf()
plt.figure(figsize=(7,5))
xt,yt=returnvals(actdata1)
plt.plot(xt,yt,'o',label='Single node data')

x=[]
y=[]
for i in range(len(sorteddata[0])):
	x.append(sorteddata[0][i][0])
	yval=np.mean([sorteddata[0][i][1],sorteddata[1][i][1],sorteddata[2][i][1],sorteddata[3][i][1],sorteddata[4][i][1]])
	y.append(yval)

plt.plot(x,y,'o',label='8 node data')
plt.legend(loc=1)
plt.ylabel('Time (seconds)',size=18)
plt.xlabel('Threads',size=18)
plt.savefig('MultiANDSingleBIGtime.pdf')









