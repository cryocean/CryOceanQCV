import os, subprocess, shlex
import sys
import datetime
path_report = "/noc/users/cryo/QCV_Cryo2/code/gen_daily_report/"
path2 = "/noc/users/cryo/QCV_Cryo2/code/"
sys.path.extend([path_report,path2])
from genReportFunTDSApril2017 import genReportFunTDS
from ftp_upDayReport import ftp_upDayReport



# copy of genAllReports.py for processing TDS (Test Data Set) for C2 data March/April 2017
# only for IOP not FDM




#matlab_path = "/noc/packages/linux_emt64/matlab/v2013a/bin/matlab" changed 24 Jan 2017
matlab_path = "/nerc/packages/matlab/2013a/bin/matlab"

myfun_path = "/noc/users/cryo/QCV_Cryo2/code/"

#date_FDM0 = datetime.date(2017,3,21)
#date_FDMf = datetime.date(2017,3,28)
date_IOP0 = datetime.date(2012,6,3)
date_IOPf = datetime.date(2012,6,3)



#n_days_FDM = (date_FDMf - date_FDM0).days + 1
n_days_IOP = (date_IOPf - date_IOP0).days + 1
#dates_FDM = [(date_FDM0 + datetime.timedelta(n)).strftime('%Y%m%d') for n in range(n_days_FDM)]
dates_IOP = [(date_IOP0 + datetime.timedelta(n)).strftime('%Y%m%d') for n in range(n_days_IOP)]

def main():
  #  for date_FDM in dates_FDM:
  #      print "generating FDM report for " + date_FDM
  #      # process data using matlab function
  #      cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; daily_report_v6_nodisplay('%s','%s'); quit" """\
  #      % (matlab_path,myfun_path,"SIR_FDM_L2",date_FDM)
  #      cmd_FDM = shlex.split(cmd_FDM)
  #      with open(os.devnull, 'w') as devnull:
  #          dum = subprocess.call(cmd_FDM,stdout = devnull)
  #      
  #      # convert figures from eps to jpg
  #      cmd_convFig = "mogrify -density 300 -quality 100% -colorspace RGB -format jpg \
  #      /noc/users/cryo/QCV_Cryo2/code/gen_daily_report/figures/*FDM*.eps"
  #      cmd_convFig = shlex.split(cmd_convFig)
  #      with open(os.devnull, 'w') as devnull:
  #          dum = subprocess.call(cmd_convFig,stdout = devnull)
  #      
  #      # generate report
  #      genReportFun("FDM")
#
  #      # upload report to ESA server
  #      ftp_upDayReport("FDM",date_FDM)

    for date_IOP in dates_IOP:
        print "generating IOP report for " + date_IOP
        # process data using matlab function
        cmd_IOP = """%s -nodisplay -nosplash -nodesktop -r "cd %s; daily_report_v6_nodisplay_TDS2017('%s','%s'); quit" """\
        % (matlab_path,myfun_path,"SIR_IOP_2",date_IOP)
        cmd_IOP = shlex.split(cmd_IOP)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_IOP,stdout = devnull)
        
        # convert figures from eps to jpg
        cmd_convFig = "mogrify -density 300 -quality 100% -colorspace RGB -format jpg \
        /noc/users/cryo/QCV_Cryo2/code/gen_daily_report/figures/I_O_P_TDS/*IOP*.eps"
        cmd_convFig = shlex.split(cmd_convFig)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_convFig,stdout = devnull)
        
        # generate report
        genReportFunTDS("IOP")

        # upload report to ESA server
    ########    ftp_upDayReport("IOP",date_IOP)
    

if __name__ == "__main__":
    main()
