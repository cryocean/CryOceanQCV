"""
Download tide gauge data from UHSLC
Old files are deleted and new ones are donwloaded
saves the data in a file in /Users/cryo/DATA/TG_data_CryOceanQCV/tg_data_hourly/

"""

import urllib
import gzip
import re
import os
import scipy as sp

path_in = "/Users/cryo/DATA/TG_data_CryOceanQCV/tg_data_hourly/" #path where data are saved
path2Data = 'http://uhslc-static.soest.hawaii.edu/nc/fdh/'

def getSlData():
    fln = [f for f in os.listdir(path_in) if os.path.isfile(path_in + f)] #files in our dir
    stationId = [f[0:13] for f in fln]
    flnDate = [f[13:21] for f in fln] #date of the last observation for each tide gauge

    #read file names on the server and check for updates
    try:
        furl = urllib.urlopen(path2Data)
        flnUrl = furl.read()
        flnUrl = re.findall('(OS.*nc)\"',flnUrl)
        idUrl = [f[0:13] for f in flnUrl] #Id of files on UHSLC: 'OS_UH-FDH000_'
        urlDate = [f[13:21] for f in flnUrl] #dates of files on UHSLC sever: '20150401'
        furl.close()
    except IOError:
        print "connection cannot be made"
    
    #download data if updated
    for i,x in enumerate(stationId): #loop over files in path_in, check for updates, and get Data 
        if x in idUrl:
            iurl = idUrl.index(x)
            fileId = urllib.URLopener()
            URL = path2Data + flnUrl[iurl]
            try:
                fileId.retrieve(URL,path_in + flnUrl[iurl])
                fileId.close()
                if urlDate[iurl] != flnDate[i]: #if dates do not coincide we delete the old file
                    os.remove(path_in + fln[i]) #delete old file
            except IOError:
                print "could not retrieve data"
                fileId.close()
                
if __name__ == "__main__":
    getSlData()
