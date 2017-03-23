"""
Download oscar data from PODAAC
ftp://podaac-ftp.jpl.nasa.gov/allData/oscar/preview/L4/oscar_third_deg/
Old files are deleted and new ones are downloaded
saves the data in a file in /scratch/general/cryosat/validation_data_2/ changed to 
/noc/mpoc/cryo/cryosat/validation_data_2/n in Jan 2017

"""

from ftplib import FTP
import ftplib
import gzip
import re
import os
import scipy as sp

path_in = "/noc/mpoc/cryo/cryosat/validation_data_2/" #path where data are saved
#path_in = "/scratch/general/cryosat/validation_data_2/" #path where data are saved

server = 'podaac-ftp.jpl.nasa.gov'
path2Data = '/allData/oscar/preview/L4/oscar_third_deg/'

def getOscarData(yr):
    fln = 'oscar_vel'+yr+'.nc.gz'
    print 'Retrieving Oscar file'+fln

    URL = path2Data + fln
    print URL
    lfn = path_in + fln
    lf  = open(lfn,'wb')

    try:
      f=FTP(server,'anonymous','cryo@noc.ac.uk')
      f.cwd(path2Data)
      try:
        f.retrbinary('RETR ' + fln, lf.write)
        lf.close()
      except IOError:
       print '***Error retrieving file '+lfn+' - check for disk full!'
       lf.close()
    #Catch ftp connection exceptions
    except ftplib.error_perm, Argument:
      print '*Permanent error', Argument, 'making connection to ftp server',
      server, 'as anonymous - check info'
    except ftplib.error_temp, Argument:
      print '*Permanent error', Argument, 'making connection to ftp server',
      server, 'as anonymous - try later'
    except ftplib.error_reply, Argument:
      print '*Unexpected reply', Argument, 'making connection to ftp server',
      server, 'as anonymous - try later'
    except ftplib.error_proto, Argument:
      print '*Reply not within protocol', Argument, 'making connection to ftp server',
      server, 'as anonymous- try later'
    except ftplib.all_errors, Argument:
      print '*Error', Argument, 'making connection to ftp server',
      server, 'as anonymous- try later'

    print 'Retrieve complete- unzipping'
    # uncompress file, save file (as NetCDF) and delete gzipped file
    fid = gzip.open(lfn,'r')
    fOut = open(lfn[:-3], 'wb')
    fOut.write(fid.read())
    fid.close()
    fOut.close()
    os.remove(lfn) #delete gzipped file
    print 'unzipping complete'

if __name__ == "__main__":
    getOscarData()
