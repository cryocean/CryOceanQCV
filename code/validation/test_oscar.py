import os, subprocess, shlex
import sys
import datetime
import re
import zipfile
path2Modules = "/noc/users/cryo/QCV_Cryo2/code/validation/"
sys.path.append(path2Modules)
#from get_sldata_ascii import getSlData #tide gauge data
#from get_cwsp_data import getCwspData #continuous wind data from NDBC
#from get_swh_wsp_data import getSwhWspData #SWH and wind data (less frequent than cwind) from NDBC
#from get_EN4_data import getEN4Data #T/S data from EN4 dataset
from get_oscar import getOscarData #Oscar current data from PODAAC
#matlab_path = "/Applications/MATLAB_R2012a.app/bin/matlab" #Kiko's mac
#matlab_path = "/noc/packages/linux_emt64/matlab/v2013a/bin/matlab" #Unix machine
matlab_path = "/nerc/packages/matlab/2013a/bin/matlab" #changed 24 Jan 2017

#dateToday = datetime.date.today()
#dateToday = dateToday.strftime('%Y%m%d')

def main():
    print 'running update...'  
    #getSlData()
    #getSwhWspData()
    #getCwspData()

    #cmd_ww3 = """%s -nodisplay -nosplash -nodesktop -r "cd %s; getWW3Data; quit" """\
    #% (matlab_path,path2Modules)
    #cmd_ww3 = shlex.split(cmd_ww3)
    #with open(os.devnull, 'w') as devnull:
    #    dum = subprocess.call(cmd_ww3,stdout = devnull)

    #getEN4Data()
    
    getOscarData('2017')


if __name__ == "__main__":
    main()

