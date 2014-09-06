# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:53:29 2014

@author: carl
"""

import os
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import matplotlib.patches as pat


#from matplotlib.colors import LogNorm, Normalize
#from astropy.cosmology import WMAP7
#from astropy.cosmology import FlatLambdaCDM
#from astropy import units as u
#from astropy import constants as const
#import numpy.ma as ma
#condition = np.str.rfind(alma_clover_files[:],'.fits') < 0
#alma_clover_fits_files = np.extract(condition,alma_clover_files)
#os.listdir(os.getcwd())
#np.ndenumerate(alma_clover_files)


#What are the important directories
directory = os.getcwd()
scripts = '/APEX2014_GratingCalibration'
co_data = '/2014/mce_data/20140726'



#initialize data directorys and files
co_directory = directory+co_data
co_directory_files = os.listdir(co_directory)
co_directory_files = co_directory_files
co_file_indexs = np.arange(0,63)

#------------------------------------------------------------------------------#
#Load the Cloverleaf Data
#make and array of file names for cloverleaf, only include if it is a fits file.
#Will be an array such that 
co_files = np.array([],ndmin=1)
i = 0 
for co_index in co_file_indexs:
    CHOP = 'co_%04i.chop' %co_index
    HK   = 'co_%04i.hk' %co_index
    INT  = 'co_%04i.int' %co_index
    ISD  = 'co_%04i.isd' %co_index
    RUN  = 'co_%04i.run' %co_index
    TS   = 'co_%04i.ts' %co_index
    index_files = np.array([(CHOP,HK,INT,ISD,RUN,TS)],
                           dtype=[('CHOP','|S20'),('HK', '|S20'),('INT','|S20'),('ISD','|S20'),('RUN','|S20'),('TS','|S20')])
    if i==0:
        co_files=index_files
    else:
        co_files=np.hstack((co_files,index_files))
    i=i+1

co_data = []
i = 0
for co_index in co_files:
    data = []
    k = 0
    for file_name in co_index:
        if k == 1 :
            f = open(co_directory+'/'+file_name)
            file_name_data = f.read()
        elif k == 4 :
            f = open(co_directory+'/'+file_name)
            file_name_data = f.read()
        elif k == 5 :
            f = open(co_directory+'/'+file_name)
            file_name_data = f.read()
        else:
            file_name_data = np.loadtxt(co_directory+'/'+file_name)
        data.append(file_name_data)
        k = k + 1
    co_data.append(data)
 #convert to an array
co_data = np.array(co_data)
 
n=0
for co_file in co_data[:,2]:  
    file_name = co_files['INT'][n]      
    f=plt.figure(n+1)
    f.clf()
    f.text(0.5,0.975,file_name,
           horizontalalignment='center',verticalalignment='top')
    ax = plt.subplot(1,1,1)
    #plt.pcoolormesh(dsigtotw, cmap = 'gray')
    plt.imshow(np.log10(np.sqrt(co_file**2)), cmap = plt.cm.jet, interpolation = 'nearest',vmax=6,vmin=3)
    #    imshow(log10(sqrt((data_all[n][1])**2)), cmap = cm.jet, interpolation = 'nearest',vmax=6,vmin=0)
    # all color maps have a "reverse" just append "_r" to the color
    # cmap = cm.cool
    # cmap = cm.gist_rainbow
    # cmap = cm.gist_heat
    # cmap = cm.hot
    # cmap = cm.jet
    # cmap = cm.prism
    # cmap = cm.spectral
    #plt.matshow(dsigtotw, cmap = 'gray') # can't handle masked arrays!
    ax.set_aspect('equal')
    plt.colorbar(orientation = 'horizontal')
    plt.savefig(co_directory+'/'+file_name+'.png' )
    n = n+1
   
n=0
for co_file in co_data[:,3]:  
    file_name = co_files['ISD'][n]      
    f=plt.figure(n+1)
    f.clf()
    f.text(0.5,0.975,file_name,
           horizontalalignment='center',verticalalignment='top')
    ax = plt.subplot(1,1,1)
    #plt.pcoolormesh(dsigtotw, cmap = 'gray')
    plt.imshow(np.log10(np.sqrt(co_file**2)), cmap = plt.cm.jet, interpolation = 'nearest',vmax=6,vmin=0)
    #    imshow(log10(sqrt((data_all[n][1])**2)), cmap = cm.jet, interpolation = 'nearest',vmax=6,vmin=0)
    # all color maps have a "reverse" just append "_r" to the color
    # cmap = cm.cool
    # cmap = cm.gist_rainbow
    # cmap = cm.gist_heat
    # cmap = cm.hot
    # cmap = cm.jet
    # cmap = cm.prism
    # cmap = cm.spectral
    #plt.matshow(dsigtotw, cmap = 'gray') # can't handle masked arrays!
    ax.set_aspect('equal')
    plt.colorbar(orientation = 'horizontal')
    plt.savefig(co_directory+'/'+file_name+'.png' )
    n = n+1
   
n=0
for co_file in co_data[:,2]:  
    file_name = co_files['INT'][n]      
    f=plt.figure(n+1)
    f.clf()
    f.text(0.5,0.975,file_name,
           horizontalalignment='center',verticalalignment='top')
    ax = plt.subplot(1,1,1)
    #plt.pcoolormesh(dsigtotw, cmap = 'gray')
    plt.imshow((np.sqrt((co_file/co_data[n,3])**2)), cmap = plt.cm.jet, interpolation = 'nearest',vmax=60,vmin=0)
    #    imshow(log10(sqrt((data_all[n][1])**2)), cmap = cm.jet, interpolation = 'nearest',vmax=6,vmin=0)
    # all color maps have a "reverse" just append "_r" to the color
    # cmap = cm.cool
    # cmap = cm.gist_rainbow
    # cmap = cm.gist_heat
    # cmap = cm.hot
    # cmap = cm.jet
    # cmap = cm.prism
    # cmap = cm.spectral
    #plt.matshow(dsigtotw, cmap = 'gray') # can't handle masked arrays!
    ax.set_aspect('equal')
    plt.colorbar(orientation = 'horizontal')
    plt.savefig(co_directory+'/'+'INTdivISD'+file_name+'.png' )
    n = n+1 
 
signal = (co_data[41,2]/co_data[37,2])
 
file_name = "co_0041_div_co_0037"      
f=plt.figure(1)
f.clf()
f.text(0.5,0.975,file_name,
       horizontalalignment='center',verticalalignment='top')
ax = plt.subplot(1,1,1)
#plt.pcoolormesh(dsigtotw, cmap = 'gray')
plt.imshow(signal, cmap = plt.cm.jet, interpolation = 'nearest',vmax=1.4,vmin=0)
#    imshow(log10(sqrt((data_all[n][1])**2)), cmap = cm.jet, interpolation = 'nearest',vmax=6,vmin=0)
# all color maps have a "reverse" just append "_r" to the color
# cmap = cm.cool
# cmap = cm.gist_rainbow
# cmap = cm.gist_heat
# cmap = cm.hot
# cmap = cm.jet
# cmap = cm.prism
# cmap = cm.spectral
#plt.matshow(dsigtotw, cmap = 'gray') # can't handle masked arrays!
ax.set_aspect('equal')
plt.colorbar(orientation = 'horizontal')
plt.savefig(co_directory+'/'+file_name+'.png' )
n = n+1 

signal = (co_data[23,2]/co_data[20,2])
 
file_name = "co_0023_div_co_0020"      
f=plt.figure(1)
f.clf()
f.text(0.5,0.975,file_name,
       horizontalalignment='center',verticalalignment='top')
ax = plt.subplot(1,1,1)
plt.imshow(signal, cmap = plt.cm.jet, interpolation = 'nearest',vmax=1.4,vmin=0)
ax.set_aspect('equal')
plt.colorbar(orientation = 'horizontal')
plt.savefig(co_directory+'/'+file_name+'.png' )

signal = (co_data[21,2]/co_data[20,2])
 
file_name = "co_0021_div_co_0020"      
f=plt.figure(1)
f.clf()
f.text(0.5,0.975,file_name,
       horizontalalignment='center',verticalalignment='top')
ax = plt.subplot(1,1,1)
plt.imshow(signal, cmap = plt.cm.jet, interpolation = 'nearest',vmax=1.4,vmin=0)
ax.set_aspect('equal')
plt.colorbar(orientation = 'horizontal')
plt.savefig(co_directory+'/'+file_name+'.png' )



#----Things above was the wrong days data -------------------------------------#
#----Using data from 20140724 -------------------------------------------------#

directory = os.getcwd()
scripts = '/APEX2014_GratingCalibration'
co_data = '/20140724'



#initialize data directorys and files
co_directory = directory+co_data
co_directory_files = os.listdir(co_directory)
co_directory_files = co_directory_files
co_file_indexs = np.arange(0,69)

#------------------------------------------------------------------------------#
#Load the  Data
#file names
co_files = np.array([],ndmin=1)
i = 0 
for co_index in co_file_indexs:
    CHOP = 'co65_%04i.chop' %co_index
    HK   = 'co65_%04i.hk' %co_index
    INT  = 'co65_%04i.int' %co_index
    ISD  = 'co65_%04i.isd' %co_index
    RUN  = 'co65_%04i.run' %co_index
    TS   = 'co65_%04i.ts' %co_index
    index_files = np.array([(CHOP,HK,INT,ISD,RUN,TS)],
                           dtype=[('CHOP','|S20'),('HK', '|S20'),('INT','|S20'),('ISD','|S20'),('RUN','|S20'),('TS','|S20')])
    if i==0:
        co_files=index_files
    else:
        co_files=np.hstack((co_files,index_files))
    i=i+1
#load actual data
co_data = []
i = 0
for co_index in co_files:
    data = []
    k = 0
    for file_name in co_index:
        if k == 1 :
            f = open(co_directory+'/'+file_name)
            file_name_data = f.read()
        elif k == 4 :
            f = open(co_directory+'/'+file_name)
            file_name_data = f.read()
        elif k == 5 :
            f = open(co_directory+'/'+file_name)
            file_name_data = f.read()
        else:
            file_name_data = np.loadtxt(co_directory+'/'+file_name)
        data.append(file_name_data)
        k = k + 1
    co_data.append(data)
 #convert to an array
co_data = np.array(co_data)
np.save(directory+'co65_20140724_all',co_data)

#-----------plot stuff --------------------------------------------------------#
n=0
for co_file in co_data[:,2]:  
    file_name = co_files['INT'][n]      
    f=plt.figure(1)
    f.clf()
    f.text(0.5,0.975,file_name,
           horizontalalignment='center',verticalalignment='top')
    ax = plt.subplot(1,1,1)
    plt.imshow(np.log10(np.sqrt(co_file**2)), cmap = plt.cm.jet, interpolation = 'nearest',vmax=6,vmin=3)
    ax.set_aspect('equal')
    plt.colorbar(orientation = 'horizontal')
    plt.savefig(co_directory+'/'+file_name+'.png' )
    n = n+1
    
#CO 6-5
signal1 = co_data[41,2]/co_data[38,2] #GI 1710
signal2 = co_data[21,2]/co_data[20,2] #GI 1680
signal3 = co_data[23,2]/co_data[20,2] #GI 1720
signal4 = co_data[24,2]/co_data[20,2] #GI 1680
signal5 = co_data[33,2]/co_data[31,2] #GI 1640
signal6 = co_data[35,2]/co_data[20,2] #GI 1700

co65_divided_data = [signal1,signal2,signal3,signal4,signal5,signal6]
co65_diveded_namesGI = ['41_38_GI1710', '21_20_GI1680','23_20_GI1720','24_20_GI1680','33_31_GI1640','35_20_GI1700']   


n=0
for divided_data in co65_divided_data:  
    file_name = co65_diveded_namesGI[n]      
    f=plt.figure(1)
    f.clf()
    f.text(0.5,0.975,'co65_'+file_name,
           horizontalalignment='center',verticalalignment='top')
    ax = plt.subplot(1,1,1)
    plt.imshow(divided_data, cmap = plt.cm.jet, interpolation = 'nearest',vmax=1.2,vmin=0)
    ax.set_aspect('equal')
    plt.colorbar(orientation = 'horizontal')
    plt.savefig(co_directory+'/'+'co65_'+file_name+'.png' )
    n = n+1
    


#CO 7-6
signal7 = co_data[49,2]/co_data[48,2] #GI 650 -- I am not sure if co65_0048 is really a clear, there is not mention in the notebook
signal8 = co_data[54,2]/co_data[53,2] #GI 790
#signal8 = co_data[23,2]/co_data[20,2] #GI 1720
#signal9 = co_data[24,2]/co_data[20,2] #GI 1680
#signal10 = co_data[33,2]/co_data[31,2] #GI 1640    
   
   
co76_divided_data = [signal7,signal8]
co76_diveded_namesGI = ['49_48_GI650', '54_53_GI790']  

n=0
for divided_data in co76_divided_data:  
    file_name = co76_diveded_namesGI[n]      
    f=plt.figure(1)
    f.clf()
    f.text(0.5,0.975,'co76_'+file_name,
           horizontalalignment='center',verticalalignment='top')
    ax = plt.subplot(1,1,1)
    plt.imshow(divided_data, cmap = plt.cm.jet, interpolation = 'nearest',vmax=1.2,vmin=0)
    ax.set_aspect('equal')
    plt.colorbar(orientation = 'horizontal')
    plt.savefig(co_directory+'/'+'co76_'+file_name+'.png' )
    n = n+1
    
 f = plt.figure(1)
ax = plt.subplot(111)
ax.invert_yaxis()
for wavelengths in wavelength_steps:
    py_values = []
    for px_value in px_range:
        py_values.append(cal_py(n,px_value,wavelengths,alpha(650),coeff))
    ax.plot(py_values,px_range)
