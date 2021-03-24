import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
#==============
def calcVOLCaod(volcpath):
	yy,T1=readVOLC(volcpath+'ICI5_030N_AOD_c.txt')
	print (yy[0])
	yy,T2=readVOLC(volcpath+'ICI5_030S_AOD_c.txt')
	yy,NN=readVOLC(volcpath+'ICI5_3090N_AOD_c.txt')
	yy,SS=readVOLC(volcpath+'ICI5_3090S_AOD_c.txt')
	TT=np.zeros(len(T1))
	for ijk in range(len(T1)):TT[ijk]=np.mean([T1[ijk],T2[ijk],NN[ijk],SS[ijk]])
	return yy,TT
def readVOLC(filename):
	fields = 'year,AOD'
	datatypes = ("float,float")
	nn=np.genfromtxt( filename ,names=fields, dtype=datatypes,skip_header=0)
	year=nn['year']
	volc=nn['AOD']
	print (year[0])
	lenV=1250
	yy=np.zeros(lenV)
	vv=np.zeros(lenV)
	for i in range(lenV):
		yy[i]=800+i
		vv[i]=np.mean(volc[(i*36):((i+1)*36)])
	return yy,vv
#============

mainpath=sys.argv[1]+'/'
RaidusEarth=6371.229

lat=[-85,-75,-65,-55,-45,-35,-25,-15,-5,5,15,25,35,45,55,65,75,85]
weight=np.cos(np.radians(lat)) 
weight/=np.sum(weight)
Area=np.zeros(43)
for ijk in range(len(Area)):
	Area[ijk]=4*np.pi*(RaidusEarth+9+ijk*0.5)**2
ConvFact=np.zeros([len(lat)*len(Area)])
for ilat in range(len(lat)):
	#ConvFact[ilat*43:(ilat+1)*43]= 0.5*weight[ilat]*Area
	ConvFact[ilat*43:(ilat+1)*43]= weight[ilat]*Area
ConvFact*=1e-9
ConvFact/=150

"""
f= open(mainpath+'Volcanic/Gao/IVI2LoadingLatHeight501-2000_Oct2012.txt')
reader = csv.reader(f, delimiter= ' ', skipinitialspace=True)
row_count = sum(1 for row in reader) 
f.close()
f= open(mainpath+'Volcanic/Gao/IVI2LoadingLatHeight501-2000_Oct2012.txt')
reader = csv.reader(f, delimiter= ' ', skipinitialspace=True)
print(row_count)
Firsttime=True
count=0
for row in reader:
	if Firsttime:
		AOD=np.zeros([row_count])	
		time=np.zeros([row_count])
		Firsttime=False
	time[count]=row[0]
	Loadings=row[1:]
	for ijk	in range(len(Loadings)):
		Loadings[ijk]=float(Loadings[ijk])*ConvFact[ijk]
	AOD[count]=np.sum(Loadings)
	count+=1
fw=open(mainpath+'Volcanic/Gao/IVI2_AOD_01-2000_Oct2012.txt','w')
for ijk in range(len(time)): fw.write(str(time[ijk])+' '+str(AOD[ijk])+'\n')
fw.close()
"""
header=7
f= open(mainpath+'Volcanic/Gao/IVI2LoadingLatHeight501-2000.txt')
reader = csv.reader(f, delimiter= ' ', skipinitialspace=True)
row_count = sum(1 for row in reader) -header
f.close()
f= open(mainpath+'Volcanic/Gao/IVI2LoadingLatHeight501-2000.txt')
reader = csv.reader(f, delimiter= ' ', skipinitialspace=True)
print(row_count)
Firsttime=True
count=0
count1=0
for row in reader:
	if count1>=header:
		if Firsttime:
			AOD1=np.zeros([row_count])	
			time=np.zeros([row_count])
			Firsttime=False
		time[count]=row[0]
		Loadings=row[1:]
		for ijk	in range(len(Loadings)):
			Loadings[ijk]=float(Loadings[ijk])*ConvFact[ijk]
		AOD1[count]=np.sum(Loadings)
		count+=1
	else: count1+=1
fw=open(mainpath+'Volcanic/Gao/IVI2_AOD_01-2000.txt','w')
for ijk in range(len(time)): fw.write(str(time[ijk])+' '+str(AOD1[ijk])+'\n')
fw.close()

#=======
yCrowley,aodCrowley=calcVOLCaod(mainpath+'Volcanic/Crowley/')

##========
#Want to plot annual means
#For Gao
lenV=len(time)/12
lenV=int(lenV)
yy=np.zeros(lenV)
#aa=np.zeros(lenV)
aa1=np.zeros(lenV)
for i in range(lenV):
	yy[i]=501+i
	#aa[i]=np.mean(AOD[(i*12):((i+1)*12)])
	aa1[i]=np.mean(AOD1[(i*12):((i+1)*12)])

plt.figure()
plt.subplot(211)
#plt.plot(yy,aa,'k')
plt.plot(yy,aa1,'k')
plt.plot(yCrowley,aodCrowley,'r')
plt.xlim([850,1850])
plt.ylim([0,0.6])
plt.subplot(212)
#plt.plot(yy,aa,'k')
plt.plot(yy,aa1,'k')
plt.plot(yCrowley,aodCrowley,'r')
plt.xlim([1850,2010])
plt.ylim([0,0.2])
plt.savefig(mainpath+'VolcAOD.png')




