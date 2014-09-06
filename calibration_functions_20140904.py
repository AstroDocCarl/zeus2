#-------------------------------------------------------------------
#----------------------GRATING Calibration Functions-----------------------
#-------------------------Carl Ferkinhoff 2014.04.20-----------------------
#-------------------------------------------------------------------
# The functions below along with the constants provided#
#-------------------------------------------------------------------
import numpy as np
import scipy as sp
import os
from scipy.optimize import leastsq



#---------------------------------------------------------------------------#
#constants - this should never change
#---------------------------------------------------------------------------#
#ZEUSII values
alpha_min= 54.33                #minimum grating 
alpha_max = 74.4  	              #maximum grating 
alpha_max_index = alpha_min     #angle at the MAXIMUM index
alpha_min_index = alpha_max     #angle at the MINIMUM index
Rd = 1.8                        #Grating drive ratio, 1 rotation of the grating 
                                #drive shaft moves the grating 1.8 degrees
Rg = 256.0                      #Stepper and drive gear ratio, 256 steps of the 
                                #motor gives 1 rotation of the grating drive shaft 



#---------------------------------------------------------------------------#
#variables - these might change depending on instrument setup and user preference
#---------------------------------------------------------------------------#

alpha_adj = 0.0	  	         #can be used if the grating experienced a shift in 
                                #the grating angle due to removing/reinstalling the grating. 
                                #Not really useful in ZEUS-2 as the grating index is set by 
                    			#the stepper-motor controller and that index can be programmed.
                    		  	#Ideally, once a calibration has been done, 
                                #the only thing necessary is to move the grating to put
                				#put a known line on a specific pixel, and then set the
                                #grating index in the stepper motor controller
                                # to the index expected given a calibration.


#It is recommended that after any changes to the next 4 variables, that you
#   recalculate the grating calibration coefficients. While recommended, it isn't required
#   as these just affect how the index relates to the grating angle, while the 
#   calibration fitting operates only on the grating angle and wavelength.
#   To change to have index_min = 0 (the negative hard limit) just change index_mid and
#   index_min to 0 and index_max to the index at the hard+ limit switch.
index_max = 2856              #corresponds to when the grating is at alpha_min, 
                                #   can be updated to reflect preferences in the index
                                #   
index_min = 0             #corresponds to when the grating is at alpha_max
                                #   can be updated to reflect preferences in the index

index_mid = 1480                 #This is the grating index of the approximate middle 
                                #     of the grating range. In truth, the value
                                #     here is arbiratary and can be any value in the
                                #     range of possible grating indexes.

#alph_mid determines the angle of index_mid
alpha_mid = (alpha_max_index-alpha_min_index)/(index_max-index_min)*(index_mid
                -index_min)+alpha_min_index



n_cent = 20.5   #not currently used for these specific functions below, but could be used when entering the calibration data
                #alternatively one could use the n_given_wave() function below 


#------------------------------------------------------------------------------#
#       Basica Functions that relate grating index (Ig) and angle (alpha)      #
#------------------------------------------------------------------------------#

#function which defines the grating angle alpha as a function of the grating or 
#   stepper motor index Ig
def alpha(Ig):
    alpha = alpha_mid+alpha_adj-(Ig-index_mid)*Rd/Rg
    return alpha
    
def index(alpha):
    index = (alpha_mid+alpha_adj-alpha)*Rg/Rd+index_mid
    return index

#function which defines the order of the grating given the wavelength being observed
#   additional logic will need to be added here to deal with the shortwavelength array

def n_given_wave(wavelength):
    if wavelength < 395:
        return 5
    if wavelength > 395:
        return 4
        
#------------------------------------------------------------------------------#
#       Basica Functions that relate the pixel positions (py,px) to y          #
#------------------------------------------------------------------------------#
# y is the effective pixel position along the effective axis. However, in the ZEUS2
# focal plane, we want to no the wavelength at the spectral (py) and spatial position
# here we define an additional function that realtes py and px to y
#def y(py,px):
#    y = np.sqrt(A*px**2+B*px*py+py**2)
        
#------------------------------------------------------------------------------#      
#           Functions for Calibrating the Grating                              #
#------------------------------------------------------------------------------#
#The Functions below load the calibration data, perform the fit to the data to 
#    determine the calibration constants then save the coefficients for future 
#    use, and lastly load previously saved coeeficients

#Load the calibration data that was previously are recorded
#   The functions below expect tab deliminated with the following structure
#   column 0 = grating index, these columns may not always be accurate, if
#                   the index_min, index_max, and index_mid are changed
#   column 1 = wavelength in microns
#   column 2 = spectral pixel position of the line (py)
#   column 3 = spatial pixel position of the line (px)
#   column 4 = grating angle, calculated at the time the data were taken
#               these should always be accurate
#   column 5 = grating order
def load_cal_data(data_directory,data_file_name):
    global cal_data
    fname = data_directory +'/'+ data_file_name
    cal_data = np.loadtxt(fname)

#Function will determine 6 coefficients that will calibrate the grating
#   i.e. it will allow us to relate a given grating index and pixel to a specific
#   wavlength. This function should work by all spatial positions. I.E. it will
#   simultaneously fit for both the spectral and spatial position.
#   The function assignes the coefficients to the "global" varialbe "coeff"
def cal_fit(alpha,py,px,n,wavelength):
    global coeff
    # Target function
    fitfunc = lambda c, alpha, py, px, n: 5/n*(c[0]*np.sin(np.pi/180*alpha)*(py+(c[6]*px**2+c[7]*px+c[8]))
                                            +c[1]*(py+(c[6]*px**2+c[7]*px+c[8]))**2+c[2]*(py+(c[6]*px**2
                                            +c[7]*px+c[8]))+c[3]+c[4]*np.sin(np.pi/180*alpha)
                                            +c[5]*(np.sin(np.pi/180*alpha))**2) 
    # Distance to the target function
    errfunc = lambda c, alpha, py, px, n, wavelength: fitfunc(c,alpha,py,px,n) - wavelength
    #initial guess for coefficients
    c0 = [  1.58701149e+00,  -5.36140009e-04,  -1.01251938e+00,-2.09972694e+02,   
          1.70513590e+03,  -1.24146904e+03,  0.000, -1.000, 0.000] 
    #perform a least squares fit to the calibration data and return the values
    # to the gloval variable "coeff"
    coeff, success = leastsq(errfunc, c0[:], args=(alpha,py,px,n,wavelength), maxfev=10**6)

#save the calibration        
def save_cal_coeff(save_directory,coeff_file_name,coeff):
        coeff_fname = save_directory + '/' + coeff_file_name
        np.savetxt(coeff_fname,coeff,fmt='%f',delimiter='\t') #saves coefficent values in a column

#load a previously save calibration
def load_coeff(coeff_directory,coeff_file_name):
    global coeff
    fname = coeff_directory +'/'+ coeff_file_name
    coeff = np.loadtxt(fname)
    
#------------------------------EXAMPLE-----------------------------------------#

#definte the directory of the data, calibration data file name, 
#   and the coefficient file name    
directory = os.getcwd()
file_name='co_cal_data_201407.csv'
coeff_file_name='cal_coeff_201407.txt'

#load the calibration data
load_cal_data(directory,file_name)

#Assign the columns of the calibration data to their values
alpha_g = cal_data[:,4] #column 4 is the grating angle that was determined when
                        #   the calibration data was taken
                        #   using this column allows one to reuse the calibration
                        #   data
py  = cal_data[:,2] #column 2 is the pixel position of the line
px  = cal_data[:,3] #column 2 is the pixel position of the line
n  = cal_data[:,5] #column 4 is the grating order of the calibration line
wavelength = cal_data[:,1] #column 1 is the wavelength of the calibration line

#Run a fit to the calibration data and print the resulting coefficients
cal_fit(alpha_g,py,px,n,wavelength)
print coeff

#save the calibration
save_cal_coeff(directory,coeff_file_name,coeff)

#load the just save calibration coefficients
load_coeff(directory,coeff_file_name)



#------------------------------------------------------------------------------#      
#       Functions for calculating wavelength and index given a calibration     #
#------------------------------------------------------------------------------#
#will calculate a wavelength at a given pixel, y, when given an order n,
#grating angle alpha, and calibration coefficient array
def cal_wavelength(n,alpha,py,px,coeff):
    n = float(n)
    py = float(py)
    px = float(px)
    alpha =float(alpha)
    wavelength =  5/n*(coeff[0]*np.sin(np.pi/180*alpha)*(py+(coeff[6]*px**2+coeff[7]*px+coeff[8]))
                    +coeff[1]*(py+(coeff[6]*px**2+coeff[7]*px+coeff[8]))**2
                    +coeff[2]*(py+(coeff[6]*px**2+coeff[7]*px+coeff[8]))
                    +coeff[3]
                    +coeff[4]*np.sin(np.pi/180*alpha)
                    +coeff[5]*(np.sin(np.pi/180*alpha))**2) 
    return wavelength    

#will calculate a grating index to place specific wavelength in grating order n
# on a pixel y, given a calibration coefficient array    
def cal_index(n,wavelength,py,px,coeff):
    py = float(py)
    px = float(px)
    n = float(n)
    wavelength = float(wavelength)
    a= coeff[5]
    b= coeff[0]*(py+(coeff[6]*px**2+coeff[7]*px+coeff[8]))+coeff[4]
    c= (coeff[1]*(py+(coeff[6]*px**2+coeff[7]*px+coeff[8]))**2
        +coeff[2]*(py+(coeff[6]*px**2+coeff[7]*px+coeff[8]))
        +coeff[3])-wavelength*n/5
    alpha = np.arcsin((-b+np.sqrt(b**2-4*a*c))/(2*a))*180/np.pi
    print alpha
    cal_index = round(index(alpha))
    return cal_index
    
def cal_py(n,px,wavelength,alpha,coeff):
    n = float(n)
    wavelength = float(wavelength)
    alpha =float(alpha)
    a = coeff[1]
    b = coeff[0]*np.sin(np.pi/180*alpha) + coeff[2]
    c = -(wavelength*n/5-(coeff[3]+coeff[4]*np.sin(np.pi/180*alpha)+coeff[5]*(np.sin(np.pi/180*alpha))**2))
    py = (((-b+np.sqrt(b**2-4*a*c))/(2*a)))-(coeff[6]*px**2+coeff[7]*px+coeff[8])
    return py

def cal_px(n,py,wavelength,alpha,coeff):
    n = float(n)
    wavelength = float(wavelength)
    alpha =float(alpha)
    a = coeff[1]
    b = coeff[0]*np.sin(np.pi/180*alpha) + coeff[2]
    c = -(wavelength*n/5-(coeff[3]+coeff[4]*np.sin(np.pi/180*alpha)+coeff[5]*(np.sin(np.pi/180*alpha))**2))
    ax = coeff[6]
    bx = coeff[7]
    cx = py + coeff[8] - (((-b+np.sqrt(b**2-4*a*c))/(2*a)))
    px = (-bx-np.sqrt(bx**2-4*ax*cx))/(2*ax)
    return round(px,1)
    
    
#--------------------------------EXAMPLE-------------------------------------#

cal_wavelength(4,alpha(132),7,coeff)

cal_index(4,439.46,7,coeff)