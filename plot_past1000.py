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
#==============
def calcVOLCaod(volcpath):
	yy,T1=readVOLC(volcpath+'ICI5_030N_AOD_c.txt')
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
	lenV=1250
	yy=np.zeros(lenV)
	vv=np.zeros(lenV)
	for i in range(lenV):
		yy[i]=800+i
		vv[i]=np.mean(volc[9+(i*36):9+((i+1)*36)])
	return yy,vv
#---
def calcVOLCcmip6(volcpath):
	fields = 'year,AOD'
	datatypes = ("float,float")
	filename=volcpath+'volcanic_sAOD_monthly_-50001-201912_new.csv'
	nn=np.genfromtxt( filename ,names=fields, dtype=datatypes,skip_header=1,delimiter=',')
	year=nn['year']
	volc=nn['AOD']
	print(year)
	lenV=len(year)/12
	print(len(year))
	print(lenV)
	lenV=int(lenV)
	yy=np.zeros(lenV)
	vv=np.zeros(lenV)
	for i in range(lenV):
		yy[i]=-499+i
		vv[i]=np.mean(volc[3+(i*12):3+((i+1)*12)])
	return yy,vv
def readVOLCaodGAO(volcpath):
	#filename=volcpath+'IVI2_AOD_01-2000_Oct2012.txt'
	filename=volcpath+'IVI2_AOD_01-2000.txt' #Use version 1 as this is the version used by PMIP
	fields = 'year,glob'
	datatypes = ("float,float")
	nn=np.genfromtxt( filename ,names=fields, dtype=datatypes)
	year=nn['year']
	volc=nn['glob']
	print(year)
	lenV=len(year)/12
	print(len(year))
	print(lenV)
	lenV=int(lenV)
	yy=np.zeros(lenV)
	vv=np.zeros(lenV)
	for i in range(lenV):
		yy[i]=501+i
		vv[i]=np.mean(volc[3+(i*12):3+((i+1)*12)])
	return yy,vv
#---------
def savedata(modcube,filename):
	fw=open(filename,'w')
	for ijk in range(len(modcube.data)): fw.write(str(modcube.coord('season_year').points[ijk]-1)+' '+str(modcube.data[ijk])+'\n')
	fw.close()
def loaddata(filename):
	fields = 'year,data'
	datatypes = ("float,float")
	nn=np.genfromtxt( filename ,names=fields, dtype=datatypes)
	return nn['year'],nn['data']
#######################
#OPTIONS FOR ANALYSIS
#start and endyear of anomaly period
st_base=1850
end_base=1900
#filter length
smoothlen=20
rmidx=int(smoothlen/2)
global pad
pad=True
#Path where the input data is stored and output will be written too
mainpath=sys.argv[1]+'/'
lw=0.6
processMODDATA=True
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


percentiles=[1,5,15,25,35,45,55,65,75,85,95,99]
nPercent=len(percentiles)
Tper=np.zeros([ntime,nPercent])
medianREC=np.zeros([ntime])
for itime in range(ntime):
	medianREC[itime]=np.percentile(data[itime,:],50)
	for iper in range(nPercent):
		Tper[itime,iper]=np.percentile(data[itime,:],percentiles[iper])
cmap = plt.cm.get_cmap('Greys')

fig=plt.figure()
ax1=fig.add_subplot(211)
ax1.set_position([0.15,0.1,0.8,0.5])
plt.plot([0,3000],[0,0],'k',linestyle='--')
for iper in range(nPercent-1):
	cnum=float(7.5-abs(iper-5))/8.
	print(cnum)
	rgba = cmap(cnum)
	plt.plot(time[849:-rmidx],Tper[849:-rmidx,iper],color=rgba)
	plt.fill_between(time[849:-rmidx],Tper[849:-rmidx,iper],Tper[849:-rmidx,iper+1],facecolor=rgba,color=rgba)

TperREC=Tper[849:-rmidx]
yearsREC=time[849:-rmidx]
medianREC=medianREC[849:-rmidx]

#======================
"""
#HADCRUT5 OBSERVATIONS
obs=iris.load(mainpath+'HadCRUT5/HadCRUT.5.0.1.0.analysis.ensemble_series.global.monthly.nc','blended air_temperature_anomaly over land with sea_water_temperature_anomaly')[0]     
addAuxCoords(obs)  
years,obs=calc_annmean(obs) 
obs=obs[:,:-1]
years=years[:-1]
ntime=len(years)
nens=np.shape(obs.data)[0]
obsdata=np.zeros([ntime,nens])        
         
for iens in range(nens):
	anomclim=obs[iens].extract(iris.Constraint( season_year=lambda yr: st_base+1<=yr<=end_base+1 ) )
	obsdata[:,iens]=obs.data[iens,:]-np.mean(anomclim.data)
	obsdata[:,iens]=butter_bandpass_filter(obsdata[:,iens], mi=smoothlen)

obsANN=iris.load(mainpath+'HadCRUT5/HadCRUT.5.0.1.0.analysis.ensemble_series.global.annual.nc','uncertainty from area not represented in the analysis (1 sigma)')[0]
addAuxCoords(obsANN) 
ERRclim=obsANN.extract(iris.Constraint( season_year=lambda yr: st_base+1<=yr<=end_base+1 ) )
ERRclim=np.mean(ERRclim.data) 
obserrANN=butter_bandpass_filter(obsANN.data, mi=smoothlen)


percentiles=[1,5,15,25,35,45,55,65,75,85,95,99]
nPercent=len(percentiles)
Tper=np.zeros([ntime,nPercent])
medianOBS=np.zeros([ntime])
for itime in range(ntime):
	medianOBS[itime]=np.percentile(obsdata[itime,:],50)
	covERR=obserrANN[itime]
	for iper in range(nPercent):
		obsENSerr=medianOBS[itime]-np.percentile(obsdata[itime,:],percentiles[iper])
		obsCOVerr=covERR*scipy.stats.norm.ppf(percentiles[iper]/100.)
		obsANOMerr=ERRclim*scipy.stats.norm.ppf(percentiles[iper]/100.)
		if percentiles[iper]<50:
			Tper[itime,iper]=medianOBS[itime]-np.sqrt(obsENSerr**2+obsCOVerr**2+obsANOMerr**2)
		else:
			Tper[itime,iper]=medianOBS[itime]+np.sqrt(obsENSerr**2+obsCOVerr**2+obsANOMerr**2)

cmap = plt.cm.get_cmap('Greys')

for iper in range(nPercent-1):
	#cnum=float(1.5-abs(iper-1))/2.
	cnum=float(7.5-abs(iper-5))/8.
	print(iper,cnum)
	rgba = cmap(cnum)
	plt.plot(years[25:-rmidx],Tper[25:-rmidx,iper],color=rgba)
	plt.fill_between(years[25:-rmidx],Tper[25:-rmidx,iper],Tper[25:-rmidx,iper+1],facecolor=rgba,color=rgba)	

TperOBS=Tper[25:-rmidx]
yearsOBS=years[25:-rmidx]
medianOBS=medianOBS[25:-rmidx]
"""
#======================
#CESM LME
nens=13
for iens in range(nens):
	if processMODDATA:
		mod=iris.load(mainpath+'CESM_LME/ProcessedData/tas_Amon_CESM1_past1000historical_r'+str(iens+1)+'i1p1f1_085001-200512.nc')[0]     
		addAuxCoords(mod)  
		years,mod=calc_annmean(mod) 
		mod=mod.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(mod))
		anomclim=mod.extract(iris.Constraint( season_year=lambda yr: st_base+1<=yr<=end_base+1 ) )
		mod.data-=np.mean(anomclim.data)
		if iens==0:moddata=np.zeros([len(years),nens])
		moddata[:,iens]=mod.data
		savedata(mod,mainpath+'SavedModData/tas_Amon_CESM1_past1000historical_r'+str(iens+1)+'i1p1f1_085001-200512.txt')
	else:
		years,mod=loaddata(mainpath+'SavedModData/tas_Amon_CESM1_past1000historical_r'+str(iens+1)+'i1p1f1_085001-200512.txt')
		if iens==0:moddata=np.zeros([len(years),nens])
		moddata[:,iens]=mod
	#
	moddata[:,iens]=butter_bandpass_filter(moddata[:,iens], mi=smoothlen)
	plt.plot(years[rmidx:-rmidx],moddata[rmidx:-rmidx,iens],'skyblue',lw=lw)
ensmean=np.zeros([len(years)])
for iyr in range(len(years)):
	ensmean[iyr]=np.mean(moddata[iyr,:])


#======================

#CMIP5
###
if processMODDATA:
	#Calculate Drift in GISS-E2-R simulations
	control=iris.load(mainpath+'CMIP5/ProcessedControlData/tas_Amon_GISS-E2-R_piControl_r1i1p141_085001-201212.nc')[0]
	addAuxCoords(control)  
	yr,control=calc_annmean(control) 
	control=control.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(control))
	control=control.extract(iris.Constraint( year=lambda yr: 851<=yr<=1999 ) )
	yrCONT=control.coord('season_year').points-1
	#correct for drift by fitting a third order polynomial
	fit=np.polyfit(yrCONT,control.data,3)
	drift=fit[3]+fit[2]*(yrCONT)+fit[1]*(yrCONT)**2+fit[0]*(yrCONT)**3
	drift[500:]=drift[500]
##
filenames=glob.glob(mainpath+'CMIP5/ProcessedData/*.nc')
nmods=len(filenames)
firstGISS=True
firstcrow=True
for idx,filename in enumerate(filenames):
	modname=filename[filename.rfind('/')+1:-3]+'.txt'
	if processMODDATA:
		mod=iris.load(filename)[0]     
		addAuxCoords(mod)  
		years,mod=calc_annmean(mod) 
		mod=mod.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(mod))
		if 'GISS-E2-R' in filename:
			mod=mod.extract(iris.Constraint( year=lambda yr: 851<=yr<=1999 ) )
			mod.data-=drift
			years=mod.coord('season_year').points-1			
		anomclim=mod.extract(iris.Constraint( season_year=lambda yr: st_base+1<=yr<=end_base+1 ) )
		mod.data-=np.mean(anomclim.data)
		savedata(mod,mainpath+'SavedModData/'+modname)
	else:
		years,mod=loaddata(mainpath+'SavedModData/'+modname)
	CMIP5data=mod.data
	CMIP5data=butter_bandpass_filter(CMIP5data, mi=smoothlen)
	if 'GISS-E2-R' in filename:
		if firstGISS:
			GISSdata=CMIP5data
			GISSyrs=years
			firstGISS=False
		else:
			GISSdata+=CMIP5data
		plt.plot(GISSyrs[rmidx:-rmidx],CMIP5data[rmidx:-rmidx],'lightgreen',lw=lw)
	else:
		if 'CCSM4' in filename:
			plt.plot(years[rmidx:-rmidx],CMIP5data[rmidx:-rmidx],'skyblue',lw=lw)
		elif 'BCC' in filename:
			plt.plot(years[rmidx:-rmidx],CMIP5data[rmidx:-rmidx],'skyblue',lw=lw)
		else:
			plt.plot(years[rmidx:-rmidx],CMIP5data[rmidx:-rmidx],'lightgreen',lw=lw)	
#
#CMIP6
filenames=glob.glob(mainpath+'CMIP6/ProcessedData/*.nc')
nmods=len(filenames)
for idx,filename in enumerate(filenames):
	modname=filename[filename.rfind('/')+1:-3]+'.txt'
	if processMODDATA:
		mod=iris.load(filename)[0]     
		addAuxCoords(mod)  
		years,mod=calc_annmean(mod) 
		mod=mod.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,weights=iris.analysis.cartography.area_weights(mod))
		anomclim=mod.extract(iris.Constraint( season_year=lambda yr: st_base+1<=yr<=end_base+1 ) )
		mod.data-=np.mean(anomclim.data)
		savedata(mod,mainpath+'SavedModData/'+modname)
	else:
		years,mod=loaddata(mainpath+'SavedModData/'+modname)
	CMIP6data=mod.data
	CMIP6data=butter_bandpass_filter(CMIP6data, mi=smoothlen)
	plt.plot(years[rmidx:-rmidx],CMIP6data[rmidx:-rmidx],'r',lw=lw)


plt.plot(yearsREC,TperREC[:,1],'dimgrey')
plt.plot(yearsREC,medianREC,'k')
plt.plot(yearsREC,TperREC[:,10],'dimgrey')
#======================
plt.plot([-1000,-2000],[-1000,-2000],'r',label='CMIP6\n(TS17)')
plt.plot([-1000,-2000],[-1000,-2000],'lightgreen',label='CMIP5\n(CU12)')
plt.plot([-1000,-2000],[-1000,-2000],'skyblue',label='CMIP5\n(GRA08)')

for iper in range(nPercent-1): 
	cnum=float(7.5-abs(iper-5))/8. 
	print(iper,cnum) 
	rgba = cmap(cnum) 
	plt.fill_between([865,945],[0.73+0.04*iper],[0.77+0.04*iper],facecolor=rgba,color=rgba)
plt.plot([860,950],[0.72+0.04*1,0.72+0.04*1],'dimgrey')
plt.plot([860,950],[0.72+0.04*10,0.72+0.04*10],'dimgrey')
plt.plot([860,950],[0.72+0.04*5.5,0.72+0.04*5.5],'k')
#======================
plt.text(955,0.72,'0.05',fontsize=9.25) 
#plt.text(1675,0.835,'0.3',fontsize=9.25) 
plt.text(955,0.899,'0.5',fontsize=9.25) 
#plt.text(1675,0.955,'0.7',fontsize=9.25) 
plt.text(955,1.09,'0.95',fontsize=9.25) 
#plt.text(1020,0.9,'HadCRUT5\nPAGES2K',fontsize=9.25) 
plt.text(1020,0.92,'PAGES2K',fontsize=9.25) 

plt.xlim([850,2010])
plt.xticks([1000,1200,1400,1600,1800,2000],['1000','1200','1400','1600','1800','2000'])
plt.ylim([-0.6,1.2])
plt.ylabel('GMST anomaly with\nrespect to 1850-1900 ($^\circ$C)')
plt.xlabel('Years (CE)')
plt.legend(loc=[0.325,0.8],frameon=False,fontsize=9.25,ncol=3)
###############################################################
####################VOLCANIC FORCING###########################
###############################################################
ax2=fig.add_subplot(212)
ax2.set_position([0.15,0.6,0.8,0.25])
plt.plot([0,3000],[0,0],'k',linestyle='--')
ax2.tick_params(axis="x",direction="in", pad=-22)


yCMIP6,aodCMIP6=calcVOLCcmip6(mainpath+'Volcanic/CMIP6/')
forceCMIP6=0.2582+aodCMIP6*-20.0
ss=np.where(yCMIP6==st_base)
ff=np.where(yCMIP6==end_base) 
anom=np.mean(forceCMIP6[ss[0][0]:ff[0][0]+1])
forceCMIP6=butter_bandpass_filter(forceCMIP6, mi=smoothlen)
plt.plot(yCMIP6[rmidx:-rmidx],forceCMIP6[rmidx:-rmidx]-anom,'r',label='TS17')

yCrowley,aodCrowley=calcVOLCaod(mainpath+'Volcanic/Crowley/')
forceCrowley=0.2582+aodCrowley*-20.0
ss=np.where(yCrowley==st_base)
ff=np.where(yCrowley==end_base) 
anom=np.mean(forceCrowley[ss[0][0]:ff[0][0]+1])
forceCrowley=butter_bandpass_filter(forceCrowley, mi=smoothlen)
plt.plot(yCrowley[rmidx:-rmidx],forceCrowley[rmidx:-rmidx]-anom,'lightgreen',label='CU12')

yGao,aodGao=readVOLCaodGAO(mainpath+'Volcanic/Gao/')
forceGao=0.2582+aodGao*-20.0
ss=np.where(yGao==st_base)
ff=np.where(yGao==end_base) 
anom=np.mean(forceGao[ss[0][0]:ff[0][0]+1])
forceGao=butter_bandpass_filter(forceGao, mi=smoothlen)
plt.plot(yGao[rmidx:-rmidx],forceGao[rmidx:-rmidx]-anom,'skyblue',label='GRA08')

plt.legend(loc="upper left",frameon=False,fontsize=9.25,ncol=5)

plt.xlim([850,2010])
plt.xticks([1000,1200,1400,1600,1800,2000],[' ',' ',' ',' ',' ',' '])
plt.ylim([-2.1,1.1])
plt.yticks([-1.5,-1.0,-0.5,0],[' ','-1.0',' ','0.0'])
plt.ylabel('Volcanic forcing\n(Wm$^{-2}$)')
plt.savefig(mainpath+'/plot_past1000.png')
#plt.savefig(mainpath+'/plot_past1000.eps', format='eps')

