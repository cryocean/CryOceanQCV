"""
Program: ftp_syncCryo.py
Author: Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
Last date modified: 13 Aug 2014

This script keeps local folders in sync with remote FTP ESA server storing
Cryosat-2 data. It keeps the same directory tree as that found in the ESA server
by creating new directories if needed. It checks for missing files or existing local files
whose size does not coincide with that of the files in the ESA server (either because
previous downloadings failed or because the files were updated) and download them.

"""
import os, os.path
import ftplib

HOST = 'science-pds.cryosat.esa.int' #remote hostname
USER = 'cryosat353' #ftp username
PASSWD = 'NoHtpQvL' #ftp password
data_sets = ('SIR_FDM_L2/LATEST',) #tupple storing the datasets
rootR = '/' 
rootL = '/noc/mpoc/cryo/cryosat' #local path where directory tree will be created
# updated Jan 2017
#rootL = '/scratch/general/cryosat' #local path where directory tree will be created


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
        ftp.cwd(rootR + dataSet) #set remote current work directory
        print "finding missing files (may take a while)...\n"
        rDirs,missFiles = listDir(ftp)
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
	        
    ftp.close() #close ftp connection
    os.chdir(rootL) #set local current directory to rootL
	        
def listDir(ftp): #recursive function
    listd = []
    fileName = []
    rDir = []
    ftp.retrlines('LIST', listd.append) #get information for all files
    listd = [x.split() for x in listd]
    dirFlag = [(x[0],x[4],x[-1]) for x in listd]
    cwd_lag = ftp.pwd() + '/'
    if dirFlag[0][0][0] == 'd': #if directory change to it, else check for missing files
        for x in dirFlag:
            ftp.cwd(ftp.pwd() + '/' + x[2])
            temp = listDir(ftp) #recursion
            fileName.extend(temp[1])
            rDir.extend(temp[0])
            ftp.cwd('..')
    else:
        missFile = findFile(cwd_lag,[(x[2],x[1]) for x in dirFlag])
        if len(missFile) != 0:
            fileName.append(missFile)
            rDir.append(cwd_lag)
    return rDir,fileName
  
def findFile(path,rNameSize): #find missing files or files with differing size
    if not os.path.isdir(rootL + path):
        os.makedirs(rootL + path)
        lNameSize = set([])
    else:
        fname = os.listdir(rootL + path)
        lNameSize = set([(x,str(os.path.getsize(rootL + path + '/' + x))) for x in fname])
    return [x[0] for x in list(set(rNameSize) - lNameSize)]
        
main()
