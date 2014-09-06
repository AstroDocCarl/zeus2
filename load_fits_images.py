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
scripts = '/analysis_scripts'
smm ='/alma_smmj02399'
clover='/alma_cloverleaf'
alma_data_folder_smm='/alma_smm_fits_20140408'
alma_data_folder_clover='/alma_cloverleaf_fits_20140408'
analysis_outputs='/analysis_outputs'
ivison='/data_from_ivison'



#initialize data directorys and files
alma_smm_directory = directory+smm+alma_data_folder_smm
alma_smm_files = os.listdir(alma_smm_directory)
alma_clover_directory = directory+clover+alma_data_folder_clover
alma_clover_files = os.listdir(alma_clover_directory)

#------------------------------------------------------------------------------#
#Load the Cloverleaf Data
#make and array of file names for cloverleaf, only include if it is a fits file.
alma_clover_fits_names=[]
for file_name in alma_clover_files:
    if np.str.rfind(file_name,'.fits')==-1:
        #print np.str.rfind(file_name,'.fits')
        continue
    else:
        #print np.str.rfind(file_name,'.fits')
        #print file_name
        alma_clover_fits_names.append(file_name)
alma_clover_fits_names=np.sort(np.array(alma_clover_fits_names))

#make an array of all of the fits tables and data
alma_clover_fits_files=np.array([[],[]],ndmin=1)
for file_name in alma_clover_fits_names:
    fits_data = np.expand_dims(np.array(pyfits.open(alma_clover_directory+'/'+file_name)),axis=1)
    if np.shape(fits_data)[0]==1:
        fits_data = np.vstack((fits_data,np.array([[0]])))
    alma_clover_fits_files=np.hstack((alma_clover_fits_files,
                                     fits_data))
    #print file_name, shape(fits_data), np.type(fits_data)

#onlyfiles = [ f for f in os.listdir(os.getcwd()) if os.path.isfile(os.path.join(os.getcwd(),f)) ]
#
# 'Clover-cont_0.0briggs.image.fits' 
# 'Clover-cont_briggs-0.5.image.fits'
# 'Clover-cont_halfArcsecUVtapered.image.fits'
# 'Clover-cont_natural.image.fits' 
# 'Clover-cont_oneArcsecUVtapered_natural.image.fits'
# 'Clover-cont_uniform.image.fits' 
# 'Clover-cont_uniform_v2.image.fits'
#
#SIMPLE  =                    T /Standard FITS                                   
#BITPIX  =                  -32 /Floating point (32 bit)                         
#NAXIS   =                    4                                                  
#NAXIS1  =                  360                                                  
#NAXIS2  =                  360                                                  
#NAXIS3  =                    1                                                  
#NAXIS4  =                    1                                                  
#EXTEND  =                    T                                                  
#BSCALE  =   1.000000000000E+00 /PHYSICAL = PIXEL*BSCALE + BZERO                 
#BZERO   =   0.000000000000E+00                                                  
#BMAJ    =   8.608750998974E-05                                                  
#BMIN    =   6.823672602574E-05                                                  
#BPA     =   4.378725051880E+01                                                  
#BTYPE   = 'Intensity'                                                           
#OBJECT  = 'Cloverleaf'                                                          
#                                                                                
#BUNIT   = 'JY/BEAM '           /Brightness (pixel) unit                         
#EQUINOX =   2.000000000000E+03                                                  
#RADESYS = 'FK5     '                                                            
#LONPOLE =   1.800000000000E+02                                                  
#LATPOLE =   1.149539000000E+01                                                  
#PC01_01 =   1.000000000000E+00                                                  
#PC02_01 =   0.000000000000E+00                                                  
#PC03_01 =   0.000000000000E+00                                                  
#PC04_01 =   0.000000000000E+00                                                  
#PC01_02 =   0.000000000000E+00                                                  
#PC02_02 =   1.000000000000E+00                                                  
#PC03_02 =   0.000000000000E+00                                                  
#PC04_02 =   0.000000000000E+00                                                  
#PC01_03 =   0.000000000000E+00                                                  
#PC02_03 =   0.000000000000E+00                                                  
#PC03_03 =   1.000000000000E+00                                                  
#PC04_03 =   0.000000000000E+00                                                  
#PC01_04 =   0.000000000000E+00                                                  
#PC02_04 =   0.000000000000E+00                                                  
#PC03_04 =   0.000000000000E+00                                                  
#PC04_04 =   1.000000000000E+00                                                  
#CTYPE1  = 'RA---SIN'                                                            
#CRVAL1  =   2.139427100000E+02                                                  
#CDELT1  =  -1.388888888889E-05                                                  
#CRPIX1  =   1.810000000000E+02                                                  
#CUNIT1  = 'deg     '                                                            
#CTYPE2  = 'DEC--SIN'                                                            
#CRVAL2  =   1.149539000000E+01                                                  
#CDELT2  =   1.388888888889E-05                                                  
#CRPIX2  =   1.810000000000E+02                                                  
#CUNIT2  = 'deg     '                                                            
#CTYPE3  = 'FREQ    '                                                            
#CRVAL3  =   6.914243032335E+11                                                  
#CDELT3  =   6.185604625459E+09                                                  
#CRPIX3  =   1.000000000000E+00                                                  
#CUNIT3  = 'Hz      '                                                            
#CTYPE4  = 'STOKES  '                                                            
#CRVAL4  =   1.000000000000E+00                                                  
#CDELT4  =   1.000000000000E+00                                                  
#CRPIX4  =   1.000000000000E+00                                                  
#CUNIT4  = '        '                                                            
#PV2_1   =   0.000000000000E+00                                                  
#PV2_2   =   0.000000000000E+00                                                  
#RESTFRQ =   6.914602886455E+11 /Rest Frequency (Hz)                             
#SPECSYS = 'TOPOCENT'           /Spectral reference frame                        
#ALTRVAL =   1.560198797112E+04 /Alternate frequency reference value             
#ALTRPIX =   1.000000000000E+00 /Alternate frequency reference pixel             
#VELREF  =                  259 /1 LSR, 2 HEL, 3 OBS, +256 Radio                 
#COMMENT casacore non-standard usage: 4 LSD, 5 GEO, 6 SOU, 7 GAL                 
#TELESCOP= 'ALMA    '                                                            
#OBSERVER= 'crazycarl42'                                                         
#DATE-OBS= '2012-07-16T23:06:34.607999'                                          
#TIMESYS = 'UTC     '                                                            
#OBSRA   =   2.139427100000E+02    ----> 14:15:46.25034                                              
#OBSDEC  =   1.149539000000E+01     ----> 11:29:43.40399999999988                                                   
#OBSGEO-X=   2.225061873184E+06                                                  
#OBSGEO-Y=  -5.440061952280E+06                                                  
#OBSGEO-Z=  -2.481682085791E+06                                                  
#DATE    = '2014-04-08T15:10:10.820000' /Date FITS file was written              
#ORIGIN  = 'CASA 4.2.0 (release r28322)' 

# 0,0 in the plot is referenced the phase center of the observation, but it doesn't exactly match NED shich gives
#Equatorial (J2000.0) 213.9426808  11.4954489 14h15m46.243s +11d29m43.62s  2.90E-01  2.90E-01    from Evans et al. 2010 (ApJS)

sigma_contour_levels=np.array([-6.0,-4.0,-2.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34,36.0,38.0])
contour_colors = ['grey','grey','grey','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k'] 
contour_styles = ['-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-']
line_thickness = 0.25
DEClimit=[-2,2]
RAlimit = [-2.0,2.0] 
plt.close('all')
minorlocator = plt.MultipleLocator(0.5)
minorformatter = plt.NullFormatter
major_beam = 'BMAJ'
minor_beam = 'BMIN'
angle_beam = 'BPA'
xy_beam    = (-1.4,-1.4)
RA_DEC_circle = [213.942655,11.495446]
RA_cent = 'CRVAL1'
DEC_cent = 'CRVAL2'
fits_files = alma_clover_fits_files
fits_file_names = alma_clover_fits_names

#component names from Alloin et al. 1997
#2D gaussian fits to the components of natural weighted image
#A (I call this 3)
#cloverleaf_fit3_2dLog.txt
#       --- ra:   186.50 +/- 0.15 pixels
#       --- dec:  174.89 +/- 0.16 pixels

#      --- ra:   186.977 +/- 0.061 pixels
#       --- dec:  174.784 +/- 0.062 pixels

#B (I call this 2)
#cloverleaf_fit2_2dLog.txt
#       --- ra:   172.80 +/- 0.15 pixels
#       --- dec:  177.80 +/- 0.15 pixels

#       --- ra:   172.245 +/- 0.047 pixels
#       --- dec:  178.051 +/- 0.047 pixels

#C (I call this 4)
#cloverleaf_fit4_2dLog.txt
#       --- ra:   196.34 +/- 0.21 pixels
#       --- dec:  188.85 +/- 0.21 pixels

#       --- ra:   196.465 +/- 0.055 pixels
#       --- dec:  189.280 +/- 0.056 pixels

#D (I call this 1)
#cloverleaf_fit1_2d_v4Log.txt
#       --- ra:   179.53 +/- 0.21 pixels
#       --- dec:  194.73 +/- 0.21 pixels

#       --- ra:   179.661 +/- 0.092 pixels
#       --- dec:  194.969 +/- 0.093 pixels

       
#Component positions in a table
components_positions = np.array([('A',    186.977, 0.061,   174.784, 0.62 , 0.0, -20.0   ),
                                 ('B',    172.245, 0.047,   178.051, 0.047, -20.0 , 0.0   ),
                                 ('C',    196.465, 0.055,   189.280, 0.056, 10.0 , 0.0    ),
                                 ('D',    179.661, 0.092,   194.969, 0.093, 0.0, 10.0    )],
                                 dtype={'names'     :['component',  'ra',   'era',  'dec',  'edec', 'ra_offset', 'dec_offset'],
                                        'formats'   :['|S10',       'float64',  'float64',  'float64',  'float64',  'float64',  'float64']})      
       
#figure, ((ax1, ax2, ax3), (ax4, ax5 , ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
f=plt.figure(1,figsize=(4,4))
ax1=plt.subplot(223)
file_name='Clover-cont_natural.image.fits'
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
oneSigma_rms = 0.000909
ax1.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks([-1.0,0,1.0],fontsize=10)
plt.xticks([-1.0,0,1.0],fontsize=10)
ax1.xaxis.set_minor_locator(minorlocator)
ax1.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=1.5,edgecolor='c',fill=False,linestyle='dashed')
ax1.add_artist(c)
#e=pat.Ellipse(xy_beam,.1,.20,45,fill=True,facecolor='r')
ax1.add_artist(e)
ax1.plot((components_positions['ra']-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0,
         (components_positions['dec']-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0,
         lw=0,c='g',marker='+',fillstyle='none',markersize=8,markeredgewidth=1)
for component in components_positions:
    ax1.annotate(component['component'],
                 ((component['ra']-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0,
                 (component['dec']-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0),
                 xytext=(component['ra_offset'],component['dec_offset']), textcoords='offset points' )    

#ax1.xaxis.set_minor_formatter(NullFormatter)
#ax1.yaxis.set_minor_formatter(NullFormatter)

ax2=plt.subplot(222)
file_name='Clover-cont_0.0briggs.image.fits' 
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
#print shape(plot_data)
oneSigma_rms = 0.0012353
ax2.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks([-1.0,0,1.0],fontsize=10)
plt.xticks([-1.0,0,1.0],fontsize=10)
ax2.xaxis.set_minor_locator(minorlocator)
ax2.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
ax2.add_artist(e)
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=1.5,edgecolor='c',fill=False,linestyle='dashed')
ax2.add_artist(c)


ax3=plt.subplot(221)
file_name='Clover-cont_uniform.image.fits' 
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
#print shape(plot_data)
oneSigma_rms = 0.00442898
ax3.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks([-1.0,0,1.0],fontsize=10)
plt.xticks([-1.0,0,1.0],fontsize=10)
ax3.xaxis.set_minor_locator(minorlocator)
ax3.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
ax3.add_artist(e)
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=1.5,edgecolor='c',fill=False,linestyle='dashed')
ax3.add_artist(c)


ax4=plt.subplot(224)
file_name='Clover-cont_halfArcsecUVtapered.image.fits'
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
#print shape(plot_data)
oneSigma_rms = 0.00148
ax4.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks([-1.0,0,1.0],fontsize=10)
plt.xticks([-1.0,0,1.0],fontsize=10)
ax4.xaxis.set_minor_locator(minorlocator)
ax4.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
ax4.add_artist(e)
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=1.5,edgecolor='c',fill=False,linestyle='dashed')
ax4.add_artist(c)

#ax5=plt.subplot(236)
#file_name='Clover-cont_oneArcsecUVtapered_natural.image.fits'
#plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
#plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
#xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
#yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
##print shape(plot_data)
#oneSigma_rms = 0.0025664
#ax5.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
#plt.ylim(DEClimit[0],DEClimit[1])
#plt.xlim(RAlimit[1],RAlimit[0])
#plt.yticks([-1.0,0,1.0],fontsize=10)
#plt.xticks([-1.0,0,1.0],fontsize=10)
#ax5.xaxis.set_minor_locator(minorlocator)
#ax5.yaxis.set_minor_locator(minorlocator)
#e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
#ax5.add_artist(e)
#xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
#c=pat.Circle(xy_circle,radius=1.5,edgecolor='c',fill=False,linestyle='dashed')
#ax5.add_artist(c)


plt.setp(plt.subplot(221).get_xticklabels(),visible=False)
plt.setp(plt.subplot(222).get_xticklabels(),visible=False)
plt.setp(plt.subplot(222).get_yticklabels(),visible=False)
#plt.setp(plt.subplot(223).get_xticklabels(),visible=False)
#plt.setp(plt.subplot(223).get_yticklabels(),visible=False)
#plt.setp(plt.subplot(224).get_yticklabels(),visible=False)
plt.setp(plt.subplot(224).get_yticklabels(),visible=False)
#plt.setp(plt.subplot(236).get_yticklabels(),visible=False)
f.subplots_adjust(hspace=0,wspace=0)
f.text(0.5,0.05,r"$\Delta$ Dec. (arcsec)",ha='center',va='center',rotation='horizontal')
f.text(0.05,0.5, r"$\Delta$ R.A. (arcsec)",ha='center',va='center',rotation='vertical')

plt.savefig('cloverleaf_all_continuum_4_v2.png',dpi=600)




#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#Load the Smm Data
#make and array of file names for cloverleaf, only include if it is a fits file.
alma_smm_fits_names=[]
for file_name in alma_smm_files:
    if np.str.rfind(file_name,'.fits')==-1:
        #print np.str.rfind(file_name,'.fits')
        continue
    else:
        #print np.str.rfind(file_name,'.fits')
        #print file_name
        alma_smm_fits_names.append(file_name)
alma_smm_fits_names=np.sort(np.array(alma_smm_fits_names))

#make an array of all of the fits tables and data
alma_smm_fits_files=np.array([[],[]],ndmin=1)
for file_name in alma_smm_fits_names:
    fits_data = np.expand_dims(np.array(pyfits.open(alma_smm_directory+'/'+file_name)),axis=1)
    if np.shape(fits_data)[0]==1:
        fits_data = np.vstack((fits_data,np.array([[0]])))
    alma_smm_fits_files=np.hstack((alma_smm_fits_files,
                                     fits_data))
    #print file_name, shape(fits_data), np.type(fits_data)

#'SMMJ02399-cont_halfArcsecUVtapered_natural_time_range.image.fits',
#       'SMMJ02399-cont_oneArcsecUVtapered_natural_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-0.2_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-0.4_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-0.6_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-0.8_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-1.0_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-1.2_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-1.4_time_range.image.fits',
#       'SMMJ02399-cont_robust_r-1.6_time_range.image.fits',
#       'SMMJ02399-cont_robust_time_range.image.fits',
#       'SMMJ02399-cont_uniform_time_range.image.fits',
#       'SMMJ02399-continuum_natural_timerange.image.fits',

#header from 'SMMJ02399-continuum_natural_timerange.image.fits'

#SIMPLE  =                    T /Standard FITS                                   
#BITPIX  =                  -32 /Floating point (32 bit)                         
#NAXIS   =                    4                                                  
#NAXIS1  =                  360                                                  
#NAXIS2  =                  360                                                  
#NAXIS3  =                    1                                                  
#NAXIS4  =                    1                                                  
#EXTEND  =                    T                                                  
#BSCALE  =   1.000000000000E+00 /PHYSICAL = PIXEL*BSCALE + BZERO                 
#BZERO   =   0.000000000000E+00                                                  
#BMAJ    =   7.273547351360E-05                                                  
#BMIN    =   7.002370225059E-05                                                  
#BPA     =   8.239395141602E+01                                                  
#BTYPE   = 'Intensity'                                                           
#OBJECT  = 'SMMJ02399-0136'                                                      
#                                                                                
#BUNIT   = 'JY/BEAM '           /Brightness (pixel) unit                         
#EQUINOX =   2.000000000000E+03                                                  
#RADESYS = 'FK5     '                                                            
#LONPOLE =   1.800000000000E+02                                                  
#LATPOLE =  -1.599667000000E+00                                                  
#PC01_01 =   1.000000000000E+00                                                  
#PC02_01 =   0.000000000000E+00                                                  
#PC03_01 =   0.000000000000E+00                                                  
#PC04_01 =   0.000000000000E+00                                                  
#PC01_02 =   0.000000000000E+00                                                  
#PC02_02 =   1.000000000000E+00                                                  
#PC03_02 =   0.000000000000E+00                                                  
#PC04_02 =   0.000000000000E+00                                                  
#PC01_03 =   0.000000000000E+00                                                  
#PC02_03 =   0.000000000000E+00                                                  
#PC03_03 =   1.000000000000E+00                                                  
#PC04_03 =   0.000000000000E+00                                                  
#PC01_04 =   0.000000000000E+00                                                  
#PC02_04 =   0.000000000000E+00                                                  
#PC03_04 =   0.000000000000E+00                                                  
#PC04_04 =   1.000000000000E+00                                                  
#CTYPE1  = 'RA---SIN'                                                            
#CRVAL1  =   3.996612500000E+01                                                  
#CDELT1  =  -1.388888888889E-05                                                  
#CRPIX1  =   1.810000000000E+02                                                  
#CUNIT1  = 'deg     '                                                            
#CTYPE2  = 'DEC--SIN'                                                            
#CRVAL2  =  -1.599667000000E+00                                                  
#CDELT2  =   1.388888888889E-05                                                  
#CRPIX2  =   1.810000000000E+02                                                  
#CUNIT2  = 'deg     '                                                            
#CTYPE3  = 'FREQ    '                                                            
#CRVAL3  =   6.468188971181E+11                                                  
#CDELT3  =   7.143401846511E+09                                                  
#CRPIX3  =   1.000000000000E+00                                                  
#CUNIT3  = 'Hz      '                                                            
#CTYPE4  = 'STOKES  '                                                            
#CRVAL4  =   1.000000000000E+00                                                  
#CDELT4  =   1.000000000000E+00                                                  
#CRPIX4  =   1.000000000000E+00                                                  
#CUNIT4  = '        '                                                            
#PV2_1   =   0.000000000000E+00                                                  
#PV2_2   =   0.000000000000E+00                                                  
#RESTFRQ =   6.459134630000E+11 /Rest Frequency (Hz)                             
#SPECSYS = 'TOPOCENT'           /Spectral reference frame                        
#ALTRVAL =  -4.202456449253E+05 /Alternate frequency reference value             
#ALTRPIX =   1.000000000000E+00 /Alternate frequency reference pixel             
#VELREF  =                  259 /1 LSR, 2 HEL, 3 OBS, +256 Radio                 
#COMMENT casacore non-standard usage: 4 LSD, 5 GEO, 6 SOU, 7 GAL                 
#TELESCOP= 'ALMA    '                                                            
#OBSERVER= 'crazycarl42'                                                         
#DATE-OBS= '2012-08-28T07:03:24.960000'                                          
#TIMESYS = 'UTC     '                                                            
#OBSRA   =   3.996612500000E+01                                                  
#OBSDEC  =  -1.599667000000E+00                                                  
#OBSGEO-X=   2.225061873184E+06                                                  
#OBSGEO-Y=  -5.440061952280E+06                                                  
#OBSGEO-Z=  -2.481682085791E+06                                                  
#DATE    = '2014-04-04T16:36:49.038000' /Date FITS file was written              
#ORIGIN  = 'CASA 4.2.0 (release r28322)'                                         

#component names from Alloin et al. 1997
#2D gaussian fits to the components of natural weighted image
#A (I call this 3)
#cloverleaf_fit3_2dLog.txt
#       --- ra:   186.50 +/- 0.15 pixels
#       --- dec:  174.89 +/- 0.16 pixels

#      --- ra:   186.977 +/- 0.061 pixels
#       --- dec:  174.784 +/- 0.062 pixels

#B (I call this 2)
#cloverleaf_fit2_2dLog.txt
#       --- ra:   172.80 +/- 0.15 pixels
#       --- dec:  177.80 +/- 0.15 pixels

#       --- ra:   172.245 +/- 0.047 pixels
#       --- dec:  178.051 +/- 0.047 pixels

#C (I call this 4)
#cloverleaf_fit4_2dLog.txt
#       --- ra:   196.34 +/- 0.21 pixels
#       --- dec:  188.85 +/- 0.21 pixels

#       --- ra:   196.465 +/- 0.055 pixels
#       --- dec:  189.280 +/- 0.056 pixels

#D (I call this 1)
#cloverleaf_fit1_2d_v4Log.txt
#       --- ra:   179.53 +/- 0.21 pixels
#       --- dec:  194.73 +/- 0.21 pixels

#       --- ra:   179.661 +/- 0.092 pixels
#       --- dec:  194.969 +/- 0.093 pixels

# 0,0 in the plot is referenced to the Valiente position, which was the phase center of the observation
sigma_contour_levels=np.array([-6.0,-4.0,-2.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34,36.0,38.0])
contour_colors = ['grey','grey','grey','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k'] 
contour_styles = ['-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-']
line_thickness = 0.25
DEClimit=[-1.5,3.5]
RAlimit = [-1.5,3.5] 
plt.close('all')
minorlocator = plt.MultipleLocator(0.5)
minorformatter = plt.NullFormatter
major_beam = 'BMAJ'
minor_beam = 'BMIN'
angle_beam = 'BPA'
xy_beam    = (2.9,2.9)
RA_DEC_circle = [39.966413263888896,-1.5994687500000002]
RA_cent = 'CRVAL1'
DEC_cent = 'CRVAL2'
fits_files = alma_smm_fits_files
fits_file_names = alma_smm_fits_names
xticks=[3.0,2.0,1.0,0,-1.0]
yticks=[3.0,2.0,1.0,0,-1.0]

#ellipse [[02:39:51.90585, -001.35.58.0875], [3.0012arcsec, 3.0000arcsec], 90.00000000deg] coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=white,
# 39.966274375000005
#-1.5994687500000002

#Fit on SMMJ02399-continuum_natural_timerange.image.fits component L1 component
#Position ---
#       --- ra:    02:39:51.818377 +/- 0.000081 s (0.001213 arcsec)  = 39.96590990416667 degrees
#       --- dec: -001.35.58.385880 +/- 0.001169 arcsec = -1.5995516333333335
#       --- ra:   195.481 +/- 0.024 pixels
#       --- dec:  188.306 +/- 0.023 pixels

#Fit on SMMJ02399-continuum_natural_timerange.image.fits component L2SW
#Position ---
#       --- ra:    02:39:51.9474 +/- 0.0013 s (0.0201 arcsec)
#       --- dec: -001.35.58.9068 +/- 0.0194 arcsec
#       --- ra:   156.79 +/- 0.40 pixels
#       --- dec:  177.89 +/- 0.39 pixels

#L1N based on HST image
#2:39;51.840 (pixels 189,213) -----> 39.966 degrees
#-1:35:57.131 -------> -1.5992030555555556 degrees
#L2 based on HST image
#2:39:52.035 (pixels 130,191) -----------> 39.9668125 degrees
#-1:35:58.249 ----------->1.5995136111111112


#Valiente et al. 2007 position, listed in ned
# Equatorial (J2000.0)   39.966125   -1.599667  02h39m51.87s  -01d35m58.8s  1.50E+00  1.50E+00    0
#       
# From Ivison et al. 2010
#L1
#L2
#L2SW
#L1N      
#Component positions in a table
components_positions = np.array([('L1',    195.481, 0.024,  188.306,    0.023,  0.0,  -20.0   ),
                                 ('L1N',   189.0,   1.0,    213.0,      1.0,  -20.0,    5.0   ),
                                 ('L2',    130.0,   1.0,    191.0,      1.0,  -15.0,    5.0    ),
                                 ('L2SW',  156.79,  0.40,   177.89,     0.39, -20.0,  -20.0    )],
                                 dtype={'names'     :['component',  'ra',   'era',  'dec',  'edec', 'ra_offset', 'dec_offset'],
                                        'formats'   :['|S10',       'float64',  'float64',  'float64',  'float64',  'float64',  'float64']})      
       
#figure, ((ax1, ax2, ax3), (ax4, ax5 , ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
f=plt.figure(1,figsize=(4,4))
ax1=plt.subplot(223)
file_name='SMMJ02399-continuum_natural_timerange.image.fits'
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
oneSigma_rms = 0.001097865
ax1.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks(yticks,fontsize=10)
plt.xticks(xticks,fontsize=10)
ax1.xaxis.set_minor_locator(minorlocator)
ax1.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=2.5,edgecolor='c',fill=False,linestyle='dashed')
ax1.add_artist(c)
#e=pat.Ellipse(xy_beam,.1,.20,45,fill=True,facecolor='r')
ax1.add_artist(e)
ax1.plot((components_positions['ra']-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0,
         (components_positions['dec']-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0,
         lw=0,c='g',marker='+',fillstyle='none',markersize=8,markeredgewidth=1)
for component in components_positions:
    ax1.annotate(component['component'],
                 ((component['ra']-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0,
                 (component['dec']-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0),
                 xytext=(component['ra_offset'],component['dec_offset']), textcoords='offset points' )    

#ax1.xaxis.set_minor_formatter(NullFormatter)
#ax1.yaxis.set_minor_formatter(NullFormatter)

ax2=plt.subplot(222)
file_name='SMMJ02399-cont_robust_time_range.image.fits' 
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
#print shape(plot_data)
oneSigma_rms = 0.001535489
ax2.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks(yticks,fontsize=10)
plt.xticks(xticks,fontsize=10)
ax2.xaxis.set_minor_locator(minorlocator)
ax2.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
ax2.add_artist(e)
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=2.5,edgecolor='c',fill=False,linestyle='dashed')
ax2.add_artist(c)
ax2.plot((components_positions['ra']-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0,
         (components_positions['dec']-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0,
         lw=0,c='g',marker='+',fillstyle='none',markersize=8,markeredgewidth=1)


ax3=plt.subplot(221)
file_name='SMMJ02399-cont_uniform_time_range.image.fits' 
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
#print shape(plot_data)
oneSigma_rms = 0.0094546
ax3.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks(yticks,fontsize=10)
plt.xticks(xticks,fontsize=10)
ax3.xaxis.set_minor_locator(minorlocator)
ax3.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
ax3.add_artist(e)
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=2.5,edgecolor='c',fill=False,linestyle='dashed')
ax3.add_artist(c)
ax3.plot((components_positions['ra']-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0,
         (components_positions['dec']-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0,
         lw=0,c='g',marker='+',fillstyle='none',markersize=8,markeredgewidth=1)


ax4=plt.subplot(224)
file_name='SMMJ02399-cont_halfArcsecUVtapered_natural_time_range.image.fits'
plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
#print shape(plot_data)
oneSigma_rms = 0.001531906
ax4.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
plt.ylim(DEClimit[0],DEClimit[1])
plt.xlim(RAlimit[1],RAlimit[0])
plt.yticks(yticks,fontsize=10)
plt.xticks(xticks,fontsize=10)
ax4.xaxis.set_minor_locator(minorlocator)
ax4.yaxis.set_minor_locator(minorlocator)
e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
ax4.add_artist(e)
xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
c=pat.Circle(xy_circle,radius=2.5,edgecolor='c',fill=False,linestyle='dashed')
ax4.add_artist(c)
#ax4.plot((components_positions['ra']-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0,
#         (components_positions['dec']-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0,
#         lw=0,c='g',marker='+',fillstyle='none',markersize=8,markeredgewidth=1)

#ax5=plt.subplot(236)
#file_name='SMMJ02399-cont_oneArcsecUVtapered_natural_time_range.image.fits'
#plot_data=fits_files[0,np.where(fits_file_names==file_name)[0][0]].data[0,0,:,:]
#plot_header=fits_files[0,np.where(fits_file_names==file_name)[0][0]].header
#xvalues=(np.arange(0,plot_header['NAXIS1'])-plot_header['CRPIX1'])*plot_header['CDELT1']*3600.0
#yvalues=(np.arange(0,plot_header['NAXIS2'])-plot_header['CRPIX2'])*plot_header['CDELT2']*3600.0
##print shape(plot_data)
#oneSigma_rms = 0.002893937
#ax5.contour(xvalues,yvalues,plot_data,levels=sigma_contour_levels*oneSigma_rms,colors=contour_colors,linewidths=line_thickness,linestyles=contour_styles)
#plt.ylim(DEClimit[0],DEClimit[1])
#plt.xlim(RAlimit[1],RAlimit[0])
#plt.yticks(yticks,fontsize=10)
#plt.xticks(xticks,fontsize=10)
#ax5.xaxis.set_minor_locator(minorlocator)
#ax5.yaxis.set_minor_locator(minorlocator)
#e=pat.Ellipse(xy_beam,width=plot_header[minor_beam]*3600,height=plot_header[major_beam]*3600,linewidth=0,angle=-plot_header[angle_beam],fill=True,facecolor='r')
#ax5.add_artist(e)
#xy_circle = ((RA_DEC_circle[0]-plot_header[RA_cent])*3600,(RA_DEC_circle[1]-plot_header[DEC_cent])*3600)
#c=pat.Circle(xy_circle,radius=2.5,edgecolor='c',fill=False,linestyle='dashed')
#ax5.add_artist(c)


plt.setp(plt.subplot(221).get_xticklabels(),visible=False)
plt.setp(plt.subplot(222).get_xticklabels(),visible=False)
plt.setp(plt.subplot(222).get_yticklabels(),visible=False)
#plt.setp(plt.subplot(234).get_yticklabels(),visible=False)
#plt.setp(plt.subplot(235).get_xticklabels(),visible=False)
plt.setp(plt.subplot(224).get_yticklabels(),visible=False)
#plt.setp(plt.subplot(236).get_yticklabels(),visible=False)
f.subplots_adjust(hspace=0,wspace=0)
f.text(0.5,0.05,r"$\Delta$ R.A. (arcsec)",ha='center',va='center',rotation='horizontal')
f.text(0.05,0.5, r"$\Delta$ Dec. (arcsec)",ha='center',va='center',rotation='vertical')

plt.savefig('smm_all_continuum_4.png',dpi=600)