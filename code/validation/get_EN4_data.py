"""
Download T/S profiles data from Metoffice EN4 data set
"""

import gzip
import re
import os
import scipy as sp
import os, subprocess, shlex
import sys
import datetime
import re
import zipfile

path2Modules = "/noc/users/cryo/QCV_Cryo2/code/validation/"

def getEN4Data():

    path_EN4 = "/noc/mpoc/cryo/cryosat/validation_data_2/EN4_TS_profiles/" #path where data are saved, updated Jan 2017
#     path_EN4 = "/scratch/general/cryosat/validation_data_2/EN4_TS_profiles/" #path where data are saved
    
    # ---------------------------- current version ----------------------------
    path2EN4Data = 'http://www.metoffice.gov.uk/hadobs/en4/data/en4-2-0/'
    FLN = "EN.4.2.0.profiles.g10."
    FLNp = "EN.4.2.0.p.profiles.g10."
    


    # --------------------------------------------------------------------------
    
    # ---------------------------- previous version ----------------------------
    """
    path2EN4Data = 'http://www.metoffice.gov.uk/hadobs/en4/data/en4-0-2/'
    FLN = "EN.4.0.2.profiles."
    FLNp = "EN.4.0.2.p.profiles."

    path2EN4Data = 'http://www.metoffice.gov.uk/hadobs/en4/data/en4-1-1/'
    FLN = "EN.4.1.1.profiles.g10."
    FLNp = "EN.4.1.1.p.profiles.g10."
    
    """
    # --------------------------------------------------------------------------
    
    flnEN4 = os.listdir(path_EN4)
    flnEN4 = [x for x in flnEN4 if ".profiles." in x]
    yearEN4 = flnEN4[-1].split('.')
    monthEN4 = yearEN4[7][4:] # month of most recent data
    yearEN4 = yearEN4[7][0:4] # year of most recent data
    dateEN4 = datetime.date(int(yearEN4),int(monthEN4),15)
    os.chdir(path_EN4)

    fln = FLN + yearEN4 + ".zip" 
    cmd_EN4 = "wget " + path2EN4Data + fln 
    cmd_EN4 = shlex.split(cmd_EN4)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_EN4,stdout = devnull)
    #unzip downloaded file (they will have '.f.' in the filename)
    if dum == 0:
        with zipfile.ZipFile(path_EN4 + fln, "r") as z:
            z.extractall(path_EN4)
        os.remove(path_EN4 + fln) #delete zip file
    
    if yearEN4 >=10:
        yearNext = int(yearEN4) + 1
    elif yearEN4 <= 3:
        yearNext = int(yearEN4) - 1
    if yearEN4 >=10 or yearEN4 <=3:
        fln2 = FLN + str(yearNext) + ".zip"
        cmd_EN4 = "wget " + path2EN4Data + fln2 
        cmd_EN4 = shlex.split(cmd_EN4)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_EN4,stdout = devnull)
        #unzip downloaded file (they will have '.f.' in the filename)
        if dum == 0:
            with zipfile.ZipFile(path_EN4 + fln2, "r") as z:
                z.extractall(path_EN4)
            os.remove(path_EN4 + fln2) #delete zip file
    
    flnEN4 = os.listdir(path_EN4) #check files that are available after downloading
    yearEN4 = [x.split('.')[7] for x in flnEN4 if ".profiles." in x]
    #download available prelimiary files (they will have '.p.' in the filename)
    os.chdir(path_EN4)
    cnt = 1
    while True:
        datei = dateEN4 + datetime.timedelta(days=(30*cnt))
        datei = datei.strftime('%Y%m')
        cnt += 1
        if datei in yearEN4:
            continue
        flni = FLNp + datei + '.nc.gz'
        cmd_EN4 = 'wget --header="Accept-Encoding: gzip" ' + path2EN4Data + flni
        cmd_EN4 = shlex.split(cmd_EN4)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_EN4,stdout = devnull)
        if dum != 0:
            break
        cmd_gzip = "gunzip " + flni 
        cmd_gzip = shlex.split(cmd_gzip)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_gzip,stdout = devnull)
        if dum != 0:
            break
    os.chdir(path2Modules)
    # change '.p.' to '.f.' in the filenames
    flnEN4 = os.listdir(path_EN4) #check files that are available after downloading
    flnEN4 = [x for x in flnEN4 if ".p." in x]
    for flni in flnEN4:
        os.rename(path_EN4 + flni,path_EN4 + flni.replace('.p.','.f.'))
                
if __name__ == "__main__":
    getEN4Data()
