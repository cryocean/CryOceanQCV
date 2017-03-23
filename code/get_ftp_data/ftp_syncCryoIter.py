"""
Program: ftp_syncCryoIter.py
Author: Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
Last date modified: 6 Jan 2017 by Chris Banks

This script keeps local folders in sync with remote FTP ESA sever storing
Cryosat-2 data. It keeps the same directory tree as is found in the remote server
by creating new directories if they are missing in the local directory tree and 
downloads all files that are missing. This sript is the same as "ftp_syncCryo.py" but
it uses an iterative method to list the directories instead of a recursive function

"""
import os, os.path
import ftplib

HOST = 'science-pds.cryosat.esa.int' #remote hostname
USER = 'cryosat353' #ftp username
PASSWD = 'NoHtpQvL' #ftp password
data_sets = ('SIR_GOP_L2',) #tupple storing the datasets that we wish to update
rootR = '/' 
rootL = '/noc/mpoc/cryo/cryosat' #local path where directory tree will be created, updated Jan 2017

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
        print "login failed"

    for dataSet in data_sets: #loop over all data_sets
        ftp.cwd(rootR + dataSet) #set remote current work directory
        print "retrieving remote directory tree (could take a while)...\n"
        rDirs = listDir(ftp,[rootR + dataSet]) #retrieve remote directory tree
        for dirs in rDirs: #loop over all directories to find missing files
            print "entering " + dirs + " and checking for missing files\n"
	    getFiles = findFile(ftp,dirs) #find missing files in directory dirs
            if len(getFiles) == 0:
                print "no missing files in " + dirs + "\n"
	    for f in getFiles: #loop over all missing files
                print "downloading file: " + f + "\n"
	        os.chdir(rootL + dirs) #set local path
	        fObj = open(f,'wb') #open file in local folder
	        ftp.retrbinary('RETR %s' %f,fObj.write) #download file
	        ct = 0
	        while os.path.getsize(f) == 0: #retrieve file again if it has size = 0
		    ct += 1
		    ftp.retrbinary('RETR %s' %f,fObj.write)
		    if ct > 10:
                        print "file " + f + "could not be retrived\n"
		        break
	        fObj.close()
	        
    ftp.close() #close ftp connection
    os.chdir(rootL) #set local current directory to rootL
	        
  
def listDir(ftp,rDir):
    dirVisit = rDir[:]
    while len(dirVisit) != 0: 
        listd = []
        ftp.cwd(dirVisit[0])
        ftp.retrlines('LIST', listd.append)
        for x in listd:
            if x[0] == 'd':
	        rDir.append(ftp.pwd() + '/' + x[56:])
                dirVisit.append(ftp.pwd() + '/' + x[56:])
            else:
                break # all branches in the directory tree end where a file is found
        dirVisit.remove(dirVisit[0])
    return rDir
  
def findFile(ftp,path):
    ftp.cwd(path)
    listd = []
    ftp.retrlines('LIST', listd.append)
    rFiles = set([x[56:] for x in listd if x[0] != 'd'])
    if not os.path.isdir(rootL + path):
        os.makedirs(rootL + path)
        lFiles = set([])
    else:
        lFiles = set(os.listdir(rootL + path))
    return list(rFiles - lFiles)
        
main()
