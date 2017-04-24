"""
Download SWH and wind speed data for each buoy.

"""

import urllib
import gzip
import re
import os
import scipy as sp
import datetime
from sets import Set
import numpy as np

def getMissingDate():
    month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    path_in = "/noc/mpoc/cryo/cryosat/validation_data_2/buoys/" # updated Jan 2017

   # path_in = "/scratch/general/cryosat/validation_data_2/buoys/"
    fileNew = os.listdir(path_in)
    fileNew = [x for x in fileNew if (x[0] == '2' and x[5] == 'c')]
    year = [int(x[0:4]) for x in fileNew]
    year = str(max(year))
    monthNew = os.listdir(path_in + fileNew[-1])
    monthNew = [x for x in monthNew if x[0] != '.']
    monthNew = list(Set(month) - Set(monthNew))

    if not monthNew: #if monthNew is empty
        monthNew = list(month)
        year = str(int(year) + 1)
    monthi = [x for x,m in enumerate(month,start=0) if m in monthNew]
    return year,monthi
    

def getCwspData():
    month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    path_in = "/noc/mpoc/cryo/cryosat/validation_data_2/buoys/" # updated Jan 2017
    #path_in = "/scratch/general/cryosat/validation_data_2/buoys/"
    year0,monthi0 = getMissingDate()
    year = year0
    monthi = list(monthi0)
    
    path_out = path_in + year + "_cwsp/"
    if not os.path.isdir(path_out):
        os.mkdir(path_out) #create directory if it does not exist
    stId = sp.loadtxt(path_in + "stationId.txt",dtype='str') # stations id
    stId = [x.lower() for x in stId] #convert to lowercase as NDBC file names are lowercase
    
    #first we check if data are in the historical folder (each file contains the whole year)
    path2Data = 'http://www.ndbc.noaa.gov/data/historical/cwind/'
    for ID in stId:
        fln = ID + 'c' + year + '.txt.gz'
        URL3 = path2Data + fln
        fileFid = urllib.URLopener()
        try:
            fileFid.retrieve(URL3,path_out + fln)
            fileFid.close()
        except IOError:
            print "No (historical) data for station " + ID 
            fileFid.close()
            continue
            
        # uncompress file, save file as txt and delete gzipped file
        fid = gzip.open(path_out + fln,'r')
        fOut = open(path_out + ID + '.txt', 'wb')
        fOut.write(fid.read())
        fid.close()
        fOut.close()
        os.remove(path_out + fln) #delete gzipped file
        
        #split content of file into months
        fid = open(path_out + ID + '.txt','r')
        h1 = fid.readline() # read header 1
        h2 = fid.readline() # read header 2
        fid.close()
        ytmp = np.loadtxt(path_out + ID + '.txt',skiprows=2)
        os.remove(path_out + ID + '.txt') #remove yearly file
        if ytmp.ndim > 1:
            for mi in monthi:
                ki = np.nonzero(ytmp[:,1] == (mi+1))
                ki = ki[0]
                yi = ytmp[ki,:]
                if not os.path.isdir(path_out + month[mi]):
                    os.mkdir(path_out + month[mi]) #create directory
                with open(path_out + month[mi] + "/" + ID + '.txt','w') as fid:
                    fid.write(h1) #write header 1
                    fid.write(h2) #write header 2
                    np.savetxt(fid,yi,fmt='%d %02d %02d %02d %02d %d %.1f %d %.1f %d') #write data in same file
                 
         
    #check if data are in the most-recent-data folder
    year,monthi = getMissingDate() 
    if (year == year0) and (monthi == monthi0): 
        for mi in monthi:
            path2Data = 'http://www.ndbc.noaa.gov/data/cwind/' + month[mi] + '/'
            try:
                furl = urllib.urlopen(path2Data)
                flnUrl = furl.read()
                fgz = re.findall('\.txt\.gz\"',flnUrl)
                ftxt = re.findall('\.txt\"',flnUrl)
                furl.close()
            except IOError:
                print "Connection cannot be made"
                continue
        
            if (not fgz) and (not ftxt): #if there are not data on the server continue
                continue
            if not os.path.isdir(path_out + month[mi]): 
                os.mkdir(path_out + month[mi]) #create directory if it does not exist  
        
            for ID in stId:
                if fgz:
                 #   print "month is " + month[mi]                    
                 #   print "ID is " + ID
                 #   print "Year is " + year
                 #   print "i is " + str(i)
                 #   fln = ID + 'b' + year + '.txt.gz'
                 #   fln = ID + str(i) + year + '.txt.gz'
                    fln = ID + str(mi+1) + year + '.txt.gz' # changed 24 April as NBDC changed naming convention

                else:
                    fln = ID  + '.txt' # changed 24 April as NBDC changed naming convention
                 #   fln = ID + 'b' + year + '.txt'  
                URL3 = path2Data + fln
                fileFid = urllib.URLopener()
                try:
                    fileFid.retrieve(URL3,path_out + month[mi] + '/' + fln)
                    fileFid.close()
                    print "got data for station " + ID
                except IOError:
                    print "No data for station " + ID 
                    fileFid.close()
                    continue
                
                # uncompress file if necessary, save file as txt and delete gzipped file
                if fgz:
                    fid = gzip.open(path_out + month[mi] + '/' + fln,'r')
                    fOut = open(path_out + month[mi] + '/' + ID + '.txt', 'wb')
                    fOut.write(fid.read())
                    fid.close()
                    fOut.close()
                    os.remove(path_out + month[mi] + '/' + fln) #delete gzipped file  
                    
    year,monthi = getMissingDate()
    monthUpdate = list(Set(monthi0) - Set(monthi))
    for x in monthUpdate:
        print "month " + month[x] + " " + year + " has been updated"
    
                        
if __name__ == "__main__":
    getCwspData()
    #import sys
    #getSwhWspData(sys.argv[1])
       
