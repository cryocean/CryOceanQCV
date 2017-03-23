"""
Program: ftp_syncCryoLatency.py
Author: Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
Last date modified: 6 Jan 2017 by Chris Banks

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
import ftplib
from time import strptime
from calendar import timegm
#HOST = 'science-pds.cryosat.esa.int' #remote hostname (baseline C)
HOST = 'baselineb.cryosat.esa.int' #remote hostname (baseline B)
USER = 'cryosat353' #ftp username (same for both baseline B and baseline C severs)
PASSWD = 'NoHtpQvL' #ftp password (same for both baseline B and baseline C severs)
data_sets = ('SIR_IOP_L2/2014/09',) #tupple storing the datasets
product = ('SIR_FDM_L2','SIR_IOP_L2','SIR_GOP_L2')
rootR = '/' 

rootL = '/noc/mpoc/cryo/cryosat' #local path where directory tree will be created, updated Jan 2017
#rootL = '/scratch/general/cryosat' #local path where directory tree will be created

sets = [x[0:10] for x in data_sets]
if product[0] in sets:
    fObjFDM = open(rootL + '/' + product[0] + '/' + "FDM_latency.txt",'a+')
if product[1] in sets:
    fObjIOP = open(rootL + '/' + product[1] + '/' + "IOP_latency.txt",'a+')
if product[2] in sets:
    fObjGOP = open(rootL + '/' + product[2] + '/' + "GOP_latency.txt",'a+')

def main():
    print "connecting to ftp server...\n" 
    try: 
        ftp = ftplib.FTP(HOST)
    except:
        print "Could not find sever\n"
    
    print "loging into ftp server...\n"
    try:
        ftp.login(USER,PASSWD)
    except:
        print "login failed\n"

    for dataSet in data_sets: #loop over all data_sets
        try:
            ftp.cwd(rootR + dataSet) #set remote current work directory
        except:
            print "Could not set current work directory"
        print "finding missing files (may take a while)...\n"
        rDirs,missFiles,latency = listDir(ftp)
        
        #Saving latency data
        if dataSet[0:10] == product[0]:
            saveLatency(fObjFDM,latency)
        if dataSet[0:10] == product[1]:
            saveLatency(fObjIOP,latency)
        if dataSet[0:10] == product[2]:
            saveLatency(fObjGOP,latency)

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
	        fObj = open(f,'wb') #open file in local folder
                ftp.retrbinary('RETR %s' %f,fObj.write) #download file
	        fObj.close()
    if product[0] in sets:
        fObjFDM.close()
    if product[1] in sets:
        fObjIOP.close()
    if product[2] in sets:
        fObjGOP.close()	        
    ftp.close() #close ftp connection
    os.chdir(rootL) #set local current directory to rootL
	        
def listDir(ftp): #recursive function
    listd = []
    fileName = []
    rDir = []
    latency = []
    ftp.retrlines('LIST', listd.append) #get information for all files
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

def saveLatency(fObj,latency):
    dataExist = [tuple(float(x) for x in y.split()) for y in fObj]
    dataExist = [' '.join([str('%.1f' %x[0]),str('%.2f' %x[1])]) for x in dataExist]
    latency = [' '.join([str('%.1f' %x[0]),str('%.2f' %x[1])]) for x in latency]
    dataNew = list(set(latency)-set(dataExist))
    for ij in dataNew:
        fObj.write(ij+'\n')

if __name__ == "__main__":
    main()
