"""
Program: ftp_getDataLatency.py
Author: Francisco Mir Calafat (francisco.calafat@noc.ac.uk) updated by Chris Banks
Last date modified: 6 March 2017

VERSION TO RETRIEVE HISTORIC BACKLOG OF DATA AND ONLY GOP DATA (can edit later if required)


This script keeps local folders in sync with remote FTP ESA server storing
Cryosat-2 data. It keeps the same directory tree as that found in the ESA server
by creating new directories if needed. This script does not update the folder 'LATEST' as it
seems to have the same data as the present month folder. It checks for missing files or existing
local files whose size does not coincide with that of the files in the ESA server (either because
previous downloadings failed or because the files were updated) and download them. This script
also computes the latency of the data, store it in a 2D matrix (time(s),latency(h)), and save the 
matrix in a file (one per product) in the product folders ('SIR_FDM_L2','SIR_IOP_L2','SIR_GOP_L2').
The latency files are updated everytime the script is run. 

"""
import os, os.path
import ftplib, socket
import time
from time import strptime
from calendar import timegm
from datetime import datetime, timedelta
tout = 1200 #timeout in seconds
HOST = 'science-pds.cryosat.esa.int' #remote hostname (baseline C)
#HOST = 'baselineb.cryosat.esa.int' #remote hostname (baseline B)
USER = 'cryosat353' #ftp username
PASSWD = 'NoHtpQvL' #ftp password
#data_sets = ('SIR_FDM_L2/2014','SIR_IOP_L2/2014','SIR_GOP_L2/2014') #data to download
product = ('SIR_FDM_L2','SIR_IOP_L2','SIR_GOP_L2')
#product = ('SIR_GOP_L2')
rootR = '/' 
# UPDATE HERE Jan 2017 
rootL = '/noc/mpoc/cryo/cryosat/testtesttest' #local path where directory tree will be created
#rootL = '/scratch/general/cryosat' #local path where directory tree will be created

# update data only for the last 2 months (FDM,IOP) and 3 months (GOP)
# to update all files for a complete year we can pass a tuple like the following instead:
# data_sets = ('SIR_FDM_L2/2014','SIR_IOP_L2/2014','SIR_GOP_L2/2014') or to update all data
# for all years we can pass data_sets = ('SIR_FDM_L2','SIR_IOP_L2','SIR_GOP_L2')
"""d1 = datetime.now()
d2 = d1.replace(day=1) - timedelta(days=2)
d3 = d2.replace(day=1) - timedelta(days=2)
data_sets = [y+ '/' + x.strftime('%Y/%m') for y in product  for x in [d1,d2,d3]]
discard = data_sets.pop(2)
discard = data_sets.pop(4)
"""
data_sets = ('SIR_GOP_L2/2010')

print "time start: " + datetime.now().strftime('%H:%M:%S')
def main():
    print "connecting to ftp server...\n"
    n_attempts = 0
    n_attempts_max = 20
    while True: 
        n_attempts += 1
        try: 
            ftp = ftplib.FTP(HOST,timeout = tout)
            connected = True
        except (ftplib.all_errors,socket.timeout):
            print "timed out while trying to connect to sever\n"
            if n_attempts >= n_attempts_max:
                raise Exception("max retry attempts reached. Exiting ...")
            print "retrying connection to server in 60s\n"
            time.sleep(60)
            connected = False
        if connected == True:
            break
    
    print "logging into ftp server...\n"
    try:
        ftp.login(USER,PASSWD,acct = "")
    except:
        print "login failed\n"
    
    for dataSet in data_sets: #loop over all data_sets
        print 'Loop at dataset ' + dataSet
        try:
            ftp.cwd(rootR + dataSet) #set remote current work directory
        except:
            continue

        print "finding missing files in " + dataSet + " (may take a while)\n"
        rDirs,missFiles,latency = listDir(ftp)       
        
        #Saving latency data
        if dataSet[0:10] == product[0]:
            fnFDM = rootL + '/' + product[0] + '/' + "FDM_latency.txt"
            saveLatency(fnFDM,latency)
        if dataSet[0:10] == product[1]:
            fnIOP = rootL + '/' + product[1] + '/' + "IOP_latency.txt"
            saveLatency(fnIOP,latency)    
        if dataSet[0:10] == product[2]:
            fnGOP = rootL + '/' + product[2] + '/' + "GOP_latency.txt"
            saveLatency(fnGOP,latency)
        print "latency data updated\n"

        if len(missFiles) == 0:
            print "no missing files " + "\n"
            continue
	for i in range(len(rDirs)): #loop over all directories with at least one missing file
            print "entering " + rDirs[i] + "\n"
            ftp.cwd(rDirs[i]) #set the remote directory
            for f in missFiles[i]:
                print "downloading file: " + f + "\n"
                os.chdir(rootL + rDirs[i]) #set local path
                if os.path.isfile(f):
                    os.remove(f)
                n_attempts = 0
                n_attempts_max = 20
                while True:
                    n_attempts += 1
                    try:
                        with open(f,'wb') as fObj: #open file in local folder
                            ftp.retrbinary('RETR %s' %f,fObj.write) #download file
                        downloaded = True
                    except (ftplib.all_errors,socket.timeout):
                        print "timed out while downloading file\n"
                        if n_attempts >= n_attempts_max:
                            raise Exception("max retry attempts reached. Exiting ...")
                        print "retrying to download file in 60s\n"
                        time.sleep(60)
                        downloaded = False
                    if downloaded == True:
                        break
    ftp.close() #close ftp connection
    os.chdir(rootL) #set local current directory to rootL
    print "time end: " + datetime.now().strftime('%H:%M:%S')
	        
def listDir(ftp): #recursive function
    listd = []
    fileName = []
    rDir = []
    latency = []
    n_attempts = 0
    n_attempts_max = 20
    while True:
        n_attempts += 1
        try:
            ftp.retrlines('LIST', listd.append) #get information for all files
            miss_found = True
        except socket.timeout:
            print "timed out while trying to find missing files\n"
            if n_attempts >= n_attempts_max:
                raise Exception("max retry attempts reached. Exiting ...")
            print "retrying to find missing files in 60s\n"
            time.sleep(60)
            miss_found = False
        if miss_found == True:
            break 
    listd = [x.split() for x in listd]
    dirFlag = [(x[0],x[4],x[-1]) for x in listd]
    cwd_lag = ftp.pwd() + '/'
    if dirFlag[0][0][0] == 'd': #if directory change to it, else check for missing files
        for x in dirFlag:
            if 'LATEST' not in x: #we do not download files in "LATEST"
                ftp.cwd(ftp.pwd() + '/' + x[2])
                temp = listDir(ftp) #recursion
                fileName.extend(temp[1])
                rDir.extend(temp[0])
                latency.extend(temp[2])
                ftp.cwd('..')
    else:
        latency.extend(getLatency(listd))
        missFile = findFile(cwd_lag,[(x[2],x[1]) for x in dirFlag])
        if len(missFile) != 0:
            fileName.append(missFile)
            rDir.append(cwd_lag)
            print rDir
    return rDir,fileName,latency
  
def findFile(path,rNameSize): #find missing files or files with differing size
    if not os.path.isdir(rootL + path):
        os.makedirs(rootL + path)
        lNameSize = set([])
    else:
        fname = os.listdir(rootL + path)
        lNameSize = set([(x,str(os.path.getsize(rootL + path + '/' + x))) for x in fname])
    return [x[0] for x in list(set(rNameSize) - lNameSize)]

def getLatency(listd): #both the timestamp of the files and the time in the file name are UTC
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    Ds = [[str(int(x[-1][19:23])+1),x[5],x[6],x[7][0:2],x[7][3:5]] if int(x[-1][23:25])>(months.index(x[5])+1) \
         else [x[-1][19:23],x[5],x[6],x[7][0:2],x[7][3:5]] for x in listd if x[-1][-3:] == 'DBL']
    Ds = [timegm(strptime(x[0]+x[1]+x[2]+x[3]+x[4],"%Y%b%d%H%M")) for x in Ds] #from UTC to seconds 
    Dm1 = [timegm(strptime(x[-1][19:23]+x[-1][23:25]+x[-1][25:27]+x[-1][28:30]+x[-1][30:32]+x[-1][32:34], \
         "%Y%m%d%H%M%S")) for x in listd if x[-1][-3:] == 'DBL']
    Dm2 = [timegm(strptime(x[-1][35:39]+x[-1][39:41]+x[-1][41:43]+x[-1][44:46]+x[-1][46:48]+x[-1][48:50], \
         "%Y%m%d%H%M%S")) for x in listd if x[-1][-3:] == 'DBL']
    Dm = [(x+y)/2.0 for x,y in zip(Dm1,Dm2)]
    latency = [(y,(x-y)/3600.0) for x,y in zip(Ds,Dm)]
    return latency    

def saveLatency(fn,latency):
    with open(fn,'r') as fObj:
        dataExist = fObj.read().splitlines()
    dataExist = [x.split() for x in dataExist]
    latency = [[str('%.1f' %x[0]),str('%.2f' %x[1])] for x in latency]
    dataExist = dict(dataExist)
    latency = dict(latency)
    dataExist.update(latency)
    dataExist = [' '.join([x,y]) for x,y in dataExist.iteritems()]
    with open(fn,'w') as fObj: 
        for ij in dataExist:
            fObj.write(ij+'\n')

if __name__ == "__main__":
    main()
