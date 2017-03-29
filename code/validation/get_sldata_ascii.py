"""
Download tide gauge data from UHSLC
Old files are deleted and new ones are donwloaded
saves the data in a file in /scratch/general/cryosat/validation_data_2/tg_data_hourly/
changed to /noc/mpoc/cryo/cryosat/validation_data_2/tg_data_hourly/ in Jan 2017

"""

import urllib
import gzip
import re
import os
import scipy as sp

path_in = "/noc/mpoc/cryo/cryosat/validation_data_2/tg_data_hourly/" #path where data are saved
#path_in = "/scratch/general/cryosat/validation_data_2/tg_data_hourly/" #path where data are saved

path2Data = 'http://ilikai.soest.hawaii.edu/woce/'

def getSlData():
    fln = [f for f in os.listdir(path_in) if os.path.isfile(path_in + f) and f[-3:] == 'dat'] #files in our dir

    for x in fln:
        fileId = urllib.URLopener()
        URL = path2Data + x
        try:
            fileId.retrieve(URL,path_in + x)
            fileId.close()
        except IOError:
            #print "could not retrieve data"
            fileId.close()

if __name__ == "__main__":
    getSlData()
