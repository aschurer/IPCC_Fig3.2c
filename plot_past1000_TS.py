import numpy as np
import csv
import matplotlib.pyplot as plt
import glob
import iris
import iris.coord_categorisation
from scipy.signal import butter, filtfilt
import scipy.stats
import sys
############################################	
def butter_bandpass_filter(data, mi=-999, ma=-999, order=2):   
	#
	if mi==-999 and ma==-999: 
		print('WARNING filter is doing nothing!')
		data_dilt=data
	else:
		if pad:
			if mi==-999:nx   = ma
			else:nx   = mi*2
			if ma==-999:nx2  = mi   
			else:nx2  = ma   
			data = np.append(np.repeat(np.mean(data[:nx2]), nx), data)
			data = np.append(data, np.repeat( np.mean(data[len(data)-nx2:]), nx))
		if ma==-999:
			b,a= butter(order, [1./float(mi)], btype='lowpass',analog=False)
		elif mi==-999:
			b,a= butter(order, [1./float(ma)], btype='highpass',analog=False)	
		else:	
			b,a= butter(order, [1./float(ma), 1./float(mi)], btype='bandpass',analog=False)
		data_filt = filtfilt(b, a, data) 
		#
		if pad:
			data_filt = data_filt[nx:-nx]  
	return data_filt
#==============
def loaddata(filename):
	fields = 'year,data'
	datatypes = ("float,float")
	nn=np.genfromtxt( filename ,names=fields, dtype=datatypes)
	return nn['year'],nn['data']
#==============
def addAuxCoords(cube):
	"""
	Add useful aux corrds to a cube
	:param cube: cube to be modified
	:return: nada as cube modified in place
	"""
	#cube.coord('longitude').circular=True # make longitude circular
	try:
		iris.coord_categorisation.add_year(cube, 'time') # add year
		iris.coord_categorisation.add_month(cube, 'time')  # add month
		iris.coord_categorisation.add_month_number(cube, 'time')  # add month
		iris.coord_categorisation.add_season_membership(cube, 'time', 'AMJJASONDJFM', name='in_AMJJASONDJFM')
		iris.coord_categorisation.add_season_year(cube, 'time', name='season_year', seasons=['AMJJASONDJFM']) 
	except (iris.exceptions.CoordinateNotFoundError,ValueError):
		pass
	for bndCoord in ['time','longitude','latitude']:
		try:
			cube.coord(bndCoord).guess_bounds()
		except (ValueError, iris.exceptions.CoordinateNotFoundError):
			pass
#==============
def calc_annmean(data):
        stmon=4
        fimon=3
        st= np.argwhere(data.coord('month_number').points==stmon)[0][0]
        fi= np.argwhere(data.coord('month_number').points==fimon)[-1][0]
        data=data[st:fi+1]
        ann_mean=data.aggregated_by(['in_AMJJASONDJFM', 'season_year'], iris.analysis.MEAN )  
        years=ann_mean.coord('season_year').points-1
        return years,ann_mean
#######################
#OPTIONS FOR ANALYSIS
#start and endyear of anomaly period
st_base=1850
end_base=1900
#filter length
smoothlen=10
rmidx=int(smoothlen/2)
global pad
pad=True
#Path where the input data is stored and output will be written too
mainpath=sys.argv[1]+'/'
#################
#
#======================
#PAGES2K RECONSTRUCTION
ensperfile=1000 # number of ensemble members for each reconstruction
#Find reconstruction files
recon_files=glob.glob(mainpath+'PAGES2K/*.txt')
nrecons=len(recon_files) # number of reconstructions
nens=nrecons*ensperfile  # total number of ensemble members
header=1 #number of header lines
#
#loop through all reconstruction files
for ifile,filename in enumerate(recon_files):
	if ifile==0:
		f= open(filename)
		reader = csv.reader(f, delimiter=' ')
		row_count = sum(1 for row in reader) 
		f.close()  
		ntime=row_count-header
		data=np.zeros([ntime,nens])	
		time=np.zeros([ntime])
	count=0
	f= open(filename)
	reader = csv.reader(f, delimiter='	')
	for row in reader:
		if count<header:
			print('Skip line',count)
		else:
			for ijk in range(len(row)):
				if ijk==0: time[count-header]=row[ijk]
				else:data[(count-header),(ifile*ensperfile)+(ijk-1)]=row[ijk]
		count+=1

st_PAGES=1961
end_PAGES=1990
for ifile,filename in enumerate(recon_files):
	en_med=np.median(data[:,ifile*ensperfile:(ifile+1)*ensperfile],axis=1)
	offset=np.mean(en_med[st_PAGES-1:end_PAGES])
	print(offset)
	data[:,ifile*ensperfile:(ifile+1)*ensperfile]-=offset

for iens in range(nens):
	data[:,iens]+=0.37
	data[:,iens]=butter_bandpass_filter(data[:,iens], mi=smoothlen)


percentiles=[5,95]
nPercent=len(percentiles)
Tper=np.zeros([ntime,nPercent])
medianREC=np.zeros([ntime])
for itime in range(ntime):
	medianREC[itime]=np.median(data[itime,:])
	for iper in range(nPercent):
		Tper[itime,iper]=np.percentile(data[itime,:],percentiles[iper])

TperREC=Tper[849:-rmidx]
yearsREC=time[849:-rmidx]
medianREC=medianREC[849:-rmidx]


modfiles=glob.glob(mainpath+'SavedModData/*.txt')
nmods=len(modfiles)
modvals=np.zeros([1149,nmods])
for imod in range(nmods):
	years,mod=loaddata(modfiles[imod])
	ss=np.where(years==851)
	ff=np.where(years==1999) 
	yrmod=years[ss[0][0]:ff[0][0]+1]
	modvals[:,imod]=butter_bandpass_filter(mod[ss[0][0]:ff[0][0]+1], mi=smoothlen)

percentiles=[5,95]
nPercent=len(percentiles)
ntime=len(yrmod)
Tper=np.zeros([ntime,nPercent])
medianMOD=np.zeros([ntime])
for itime in range(ntime):
	medianMOD[itime]=np.mean(modvals[itime,:])
	for iper in range(nPercent):
		Tper[itime,iper]=np.percentile(modvals[itime,:],percentiles[iper])

TperMOD=Tper[rmidx:-rmidx]
yearsMOD=yrmod[rmidx:-rmidx]
medianMOD=medianMOD[rmidx:-rmidx]
#
filename=mainpath+'Observations/AR6_average_temperature.csv' 
fields = 'year,glob'
datatypes = ("float,float")
nn=np.genfromtxt( filename ,names=fields, dtype=datatypes,delimiter=',',skip_header=1)
years=np.flip(nn['year'],0)
obsdata=np.flip(nn['glob'],0)
obs=butter_bandpass_filter(obsdata, mi=smoothlen)
yearsOBS=years[rmidx:-rmidx]
obs=obs[rmidx:-rmidx]
#
plt.close()
fig=plt.figure()
ax1=fig.add_subplot(211)
ax1.set_position([0.125,0.1,0.8,0.6])
plt.plot([0,3000],[0,0],'k',linestyle='--')
#
plt.fill_between(yearsREC,TperREC[:,0],TperREC[:,1],facecolor='grey',color='grey',alpha=0.45)
plt.fill_between(yearsMOD,TperMOD[:,0],TperMOD[:,1],facecolor='salmon',color='salmon',alpha=0.45)
#
plt.plot(yearsOBS,obs,'g',label='Observed')
plt.plot(yearsREC,medianREC,'k',label='Reconstructed')
plt.plot(yearsMOD,medianMOD,'r',label='Simulated')
#======================
#
plt.fill_between([1156.2,1217],[0.9],[1],facecolor='grey',color='grey',alpha=0.45)
plt.fill_between([1511.5,1572.1],[0.9],[1],facecolor='salmon',color='salmon',alpha=0.45)

plt.xlim([850,2020])
plt.xticks([1000,1200,1400,1600,1800,2000],['1000','1200','1400','1600','1800','2000'])
#
if smoothlen<=10:plt.ylim([-0.75,1.05])
else:plt.ylim([-0.6,0.9])
plt.ylabel('Global temperature ($^\circ$C)')
plt.xlabel('Years (CE)')
plt.legend(loc='upper left',frameon=False,fontsize=9.25,ncol=3)
plt.savefig(mainpath+'past1000_TS_filter'+str(smoothlen)+'.png')
#plt.savefig(mainpath+'past1000_TS_filter'+str(smoothlen)+'.eps', format='eps')
