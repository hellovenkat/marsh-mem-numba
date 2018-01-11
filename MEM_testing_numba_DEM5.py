#gdal needs to be installed on the computer
import time
import os
import subprocess
import gdal
import rasterio
import numpy as np
from numba import jit
#import threading
start = time.time()## Starting the counter
time.strftime("%H:%M:%S")
os.chdir(r"/home/venkat/Documents/files_rci/")
outputDirectory = 'DEM5/'
#make a folder in the current directory called output if the folder doesnt already exist
if not os.path.exists(outputDirectory) : os.makedirs(outputDirectory)
#Variables
#gDEM = "NI_5m.tif"
#gDEM = "Beaufort_DEM2013_meters.tif" #elevation in meters relative to navd88
#gDEM = "Beaufor_3mDEM_4.tif"
gDEM ="Small_test_area.tif"
#mask = "F:\PIE_MEM3D_Simulation\Master_Data\pie_land_mask.shp"
#Mean sea level at 2.230 m (local datum) Fort Pulaski, GA
#MHW = 1.009 #mean high water (m) (relative to mean sea level!!)
#nadv88 2.301 relative to station
#Calculated the MSL and range based on the average in the area
MSL = -0.026 #in meters in relation to nad88
AmpM = 2.23125/2 #tidal amplitude in meters
MHW = 1.089625 #relative to nad88  https://www.ngs.noaa.gov/Tidal_Elevation/diagram.jsp?PID=CK0694&EPOCH=1983-2001
MHWcm = MHW * 100  #cm
#MLW = -1.099 #mean low water (m) relative to msl
MLW = -1.141625 #relative to NAD88
MLWcm = MLW *100  #cm
MinE = MSL- 0.10 #minimum elevation in m
MaxE = MHW + 0.3 # max elevation in m
fE = MaxE + 1  #representing the maximum elevation needed for model simulation, was calculated as the maximum vegetated elevation plus the total SLR
Trange = MHWcm - MLWcm #Tidal Range
gSLR = 100 #global SLR in cm
Amp= Trange/2  #tidal amplitude
T = 100 #time to run model in years
A = 0.317 #local rate of SLR in cm/yr (Fort Pulaski Georgia Tide Gauge)
sB = (gSLR - A*T)/T**2 #accelerating term for SLR
#p for patens a for alterniflora
aa =2760.28744894422 #biomass coefficient (changing to cm if mult by 0.0001)
ab =-3201.81817532099#biomass coefficient
ac =805.089047566297 #biomass coefficient
q = 2.8 # g/g  #proportional to the rate of sediment loading (trapping efficientcy)
kr = 0.1 # g/g  refracorty fraction of OM
#Edwards chhanged to 14e-6 g/cm3 to match PIE data reported by Cavatorta et al. (2003)
#m = 21e-6 #g/cm^3
#Calculated suspended sediment based on average in area
m= 3.35e-5 #suspended sedminent concentration g/cm3 https://pubs.usgs.gov/sir/2008/5150/pdf/sir20085150.pdf  or https://www.waterqualitydata.us/portal/
f = 704 #the frequency of semi-diurnal tides per year
k1 = 0.085 #pure organic g/cm3 packing density
k2 = 1.99 #pure inorganic g/cm3
BGTR = 1 #belowground turnover rate
RSR = 2 #root-shoot ratio
outtime = 100 #years
t_int = 1
#Dictionaries
# Threading
@jit( nogil=True, cache=True)
def PlatDepth(arra, ysr, year, f_fE, f_MinE, f_AmpM, f_Amp, f_MHWcm, f_Trange, f_aa, f_ab, f_ac, f_k2, f_q, f_m, f_f, f_k1, f_kr, f_RSR, f_BGTR):
    platList = [None] * len(arra)
    #arrayYear = np.empty(shape=(f_y_pixels, f_x_pixels))
    for idx, lin in enumerate(arra):
        platList[idx] = readFn(ysr, lin, f_fE, f_MinE, f_AmpM, f_Amp, f_MHWcm, f_Trange, f_aa, f_ab, f_ac, f_k2, f_q, f_m, f_f, f_k1, f_kr, f_RSR, f_BGTR)
    #Storing the output array in the total array
    totalList[year-1] = platList

def readFn(yslr, line, r_fE, r_MinE, r_AmpM, r_Amp, r_MHWcm, r_Trange, r_aa, r_ab, r_ac, r_k2, r_q, r_m, r_f, r_k1, r_kr, r_RSR, r_BGTR):
    List=list()
    readArray = np.empty(len(line))
    ind = 0
    for pixel in line:
        if pixel== "-32768":
                readArray[ind] = "-9999"
        elif pixel== "-3.40282e+38":
                readArray[ind] = "-9999"
        elif pixel== "-9999":
                readArray[ind] = "-9999"
        elif pixel > r_fE:
                readArray[ind] = "-9999"
        elif pixel <= r_MinE:
                ZZ = float(pixel)
                Z =  ZZ*100 #change units to centimeters
                ZZamp = r_AmpM-ZZ #meters
                Zamp  = r_Amp - Z #cm
                #D = min(1, ((MHWcm) - Z)/Trange) #Dimensionless depth in cm values between 0 and 1
                #D = max(0, D)
                #absD = max(Amp - Z, 0) # abs Depth when using elevation relative to NAVD88
                #absD = round(max(MHWcm - Z, 0), 2)  #use this when using elevation relative to MSL
                #B = max((aa*absD + ab*absD*absD + ac), 0)  #When using depth to calculate biomass
                #B = (aa*Z + ab*Z**2 + ac)
                #B = max((aa*ZZamp + ab*ZZamp*ZZamp + ac), 0)  # when using elevation to calculate biomass convert g/m2 to g/cm2
                #B = 0 #g/m2 to g/cm2
                DEDT = 0    #    #v5.4=(q*m*f*(Dreal/2) + kr*B)*(LOI/k1+(1-LOI)/k2)
                #print DEDT
                #Forcing elevation into meters
                Z = (Z + ((DEDT*10)- yslr))  #HMMM maybe we need to find out how to use the generated Z's  for future instread  of same input DEM
                #Z = (Z + (DEDT- yslr))
                ZZZ = Z/100                 #put elevation output back into meters 
                readArray[ind] = ZZZ #append the depth into the list
        elif pixel > r_MinE:
                ZZ = float(pixel)
                Z =  ZZ*100 #change units to centimeters
                ZZamp = r_AmpM-ZZ #meters
                Zamp  = r_Amp - Z #cm
                D = min(1, ((r_MHWcm) - Z)/r_Trange) #Dimensionless depth in cm values between 0 and 1
                D = max(0, D)
                #absD = max(Amp - Z, 0) # abs Depth when using elevation relative to NAVD88
                absD = round(max(r_MHWcm - Z, 0), 2)  #use this when using elevation relative to MSL
                #B = max((aa*absD + ab*absD*absD + ac), 0)  #When using depth to calculate biomass
                #B = (aa*Z + ab*Z**2 + ac)
                B = max((r_aa*Zamp + r_ab*Zamp**2 + r_ac), 0)*0.0001 #g/m2 to g/cm2  # when using elevation to calculate biomass convert g/m2 to g/cm2
                DEDT = (((1/r_k2)*(r_q*r_m*r_f*absD*0.5*D))+((1/r_k1)*(r_kr*r_RSR*r_BGTR*B)))     #    #v5.4=(q*m*f*(Dreal/2) + kr*B)*(LOI/k1+(1-LOI)/k2)
                #print DEDT
                #Forcing elevation into meters
                Z = (Z + ((DEDT*10)- yslr))  #HMMM maybe we need to find out how to use the generated Z's  for future instread  of same input DEM
                #print Z
                #Z = (Z + (DEDT- yslr))
                ZZZ = Z/100                 #put elevation output back into meters
                readArray[ind] = ZZZ #append the depth into the list
        ind=ind+1
    return readArray


#########################
#Create raster from array
#########################

def array_to_raster(array, Output_file, x_pixels, y_pixels, PIXEL_SIZE, x_min, y_max):
    dst_filename = Output_file
    #print(dst_filename)
    #Get values
    #Number of rows. In this case I'm taking the information from my input raster.
    #x_pixels = raster.width
    #print(x_pixels)
    #Number of columns. Also taking the information from my input raster.
    #y_pixels = raster.height
    #print(y_pixels)
    #Size of the pixels, pretty sure the default unit is meteres. Also taking this info from input raster.
    #PIXEL_SIZE = raster1.GetGeoTransform()[1]
    #print(PIXEL_SIZE)
    #x_min and y_max are the longitude and latitude values of the top left corner of the image. The "GetGeoTransform" tool gets longitude and latitude
    #information for each corner of the image.
    #x_min = raster1.GetGeoTransform()[0]
    #print(x_min)
    #y_max = raster1.GetGeoTransform()[3]
    #print(y_max)
    #This is the projection information. WKT stands for "Well-known text". If you input the projection manually you need the WKT name for the projection.
    #wkt_projection = raster1.GetProjection()
    
    #This is the type of output file you want. Probably always going to be geotiff.
    driver = gdal.GetDriverByName('GTiff')

    #This is setting the variable "dataset" to the driver.Create command, which creates new rasters based on your input info.
    dataset = driver.Create(
        dst_filename,
        x_pixels,
        y_pixels,
        1,
        gdal.GDT_Float32, )
    if dataset != None:
        #I think this is setting the pixel size in the output raster.
        dataset.SetGeoTransform((
            x_min,    # 0
            PIXEL_SIZE,  # 1
            0,                      # 2
            y_max,    # 3
            0,                      # 4
            -PIXEL_SIZE))

        #This sets the projection info for the output raster.
        dataset.SetProjection(wkt_projection)
        #set no data value
        NoData_value = -9999
	myarray = np.asarray(array)
        dataset.GetRasterBand(1).WriteArray(myarray)
        dataset.FlushCache()  # Write to disk.
        return dataset, dataset.GetRasterBand(1)  #If you need to return, remenber to return  also the dataset because the band doesnt live without dataset.
        #assign no data value
    dataset.SetNoDataValue(NoData_value)

#Setting sea levels
SL = list(range(0, T+1)) # converted to list object as range object does not support assignment
#SL[0]=A#Sea level at time 0 is referenced to the current 0cm, not 0.32cm above it as sea level has not risen yet. This happens in year 1
for i in range(0, T+1):
    SL[i] = A*i + sB * i**2
ySLR = list(range(0, T+1)) # converted to list object as range object does not support assignment
for i in range(1, T+1):
    ySLR[i] = round(SL[i] - SL[i - 1],5)
print ("sea levels finished")


dem_input = gDEM
#Open raster in both gdal and rasterio
ras= rasterio.open(dem_input)
ras1= gdal.Open(dem_input)
#Create empty list with same shape as input raster, file1, and fill with 0s
year1_list= np.empty(shape=(ras.height, ras.width))
x_pixels = ras.width
#Number of columns. Also taking the information from my input raster.
y_pixels = ras.height
#Intializing the global variables so that they need not be passed while writing to a file
#Size of the pixels, pretty sure the default unit is meteres. Also taking this info from input raster.
PIXEL_SIZE = ras1.GetGeoTransform()[1]
#x_min and y_max are the longitude and latitude values of the top left corner of the image. The "GetGeoTransform" tool gets longitude and latitude
#information for each corner of the image.
x_min = ras1.GetGeoTransform()[0]
y_max = ras1.GetGeoTransform()[3]
wkt_projection = ras1.GetProjection()

#Declaring the Array which contains all the arrays for every year
totalList = np.empty((T,y_pixels,x_pixels))

############
#Running the 1st year directly from the .tif file 
############
x = 0
for lin in ras.read(1):
    year1_list[x] = (readFn(ySLR[1], lin, fE, MinE, AmpM, Amp, MHWcm, Trange, aa, ab, ac, k2, q, m, f, k1, kr, RSR, BGTR))
    x = x + 1

#print(year1_list)
#Storing the result of first year('year1_list') into the totalArray
totalList[0] = (year1_list)


#Running all the years by passing the output array of the previous year as a parameter
for i in range(2, T+t_int, t_int):
    PlatDepth(totalList[i-2], ySLR[i], i, fE, MinE, AmpM, Amp, MHWcm, Trange, aa, ab, ac, k2, q, m, f, k1, kr, RSR, BGTR)

#Writing the files from total array
for i in range(1, T+t_int, t_int):
    name_i = str(i)
    outputFileName = outputDirectory + '//' + 'year_' + name_i + '.tif'
    array_to_raster(totalList[i-1], outputFileName, x_pixels, y_pixels, PIXEL_SIZE, x_min, y_max)




end = time.time()
print ("Time Elapsed (milli sec):")
print ((end-start)*1000)
print ("Time Elapsed (min):")
print ((end - start)/60)
print ("Time Elapsed (hour):")
print ((end - start)/3600)



