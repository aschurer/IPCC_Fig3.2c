"""
Script needed to correct timestamps in CESM files
So correct dates get read in by iris load commands
"""
import iris
import glob
import os
import numpy as np
import sys

inputdir=sys.argv[1]+'/CESM_LME/RawData/'
outputdir=sys.argv[1]+'/CESM_LME/ProcessedData/'
os.system('mkdir -p '+outputdir)

std_name='air_temperature'

nens=13
for iens in range(nens):
	os.system('rm -rf '+inputdir+'temp.nc')
	os.system('rm -rf '+inputdir+'temp1.nc')
	if iens <=8:ensname='0'+str(iens+1)
	else:ensname=str(iens+1)
	command='ncrcat '+inputdir+'b.e11.BLMTRC5CN.f19_g16.0'+ensname+'.cam.h0.TREFHT.085001-184912.nc '+inputdir+'b.e11.BLMTRC5CN.f19_g16.0'+ensname+'.cam.h0.TREFHT.185001-200512.nc '+inputdir+'temp.nc'
	print(command)
	os.system(command)
	cube=iris.load(inputdir+'temp.nc','Reference height temperature')[0]
	cube.standard_name = std_name
	print(cube)
	filename=inputdir+'temp1.nc'
	iris.fileformats.netcdf.save(cube, filename, netcdf_format='NETCDF4')
	outname=outputdir+'tas_Amon_CESM1_past1000historical_r'+str(iens+1)+'i1p1f1_085001-200512.nc'	
	command="ncap2 -s time-=15 "+filename+" "+outname
	print(command)
	os.system(command)
	command="rm -rf "+filename
	print(command)
	os.system(command)
os.system('rm -rf '+inputdir+'temp.nc')
os.system('rm -rf '+inputdir+'temp1.nc')
