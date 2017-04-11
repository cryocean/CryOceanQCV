import os, subprocess, shlex
import sys
import datetime
path_report = "/noc/users/cryo/QCV_Cryo2/code/gen_daily_report/"
path2 = "/noc/users/cryo/QCV_Cryo2/code/"
sys.path.extend([path_report,path2])
from genReportFun import genReportFun
from ftp_upDayReport import ftp_upDayReport

#matlab_path = "/noc/packages/linux_emt64/matlab/v2013a/bin/matlab" changed 24 Jan 2017
matlab_path = "/nerc/packages/matlab/2013a/bin/matlab"

myfun_path = "/noc/users/cryo/QCV_Cryo2/code/"

configDateFile = "/noc/users/cryo/QCV_Cryo2/code/dailyconfig.txt"  # defines start and stop dates


#date_FDM0 = datetime.date(2016,1,4)
#date_IOP0 = datetime.date(2015,12,31)
#date_FDMf = datetime.date(2016,1,4)
#date_IOPf = datetime.date(2015,12,31)

# to use config file dailyconfig.txt use split to change the arguments on a given line

# read start and stop dates from file specified as configDateFile
i=0
f=open(configDateFile)
for line in f:
    if i == 0:
        t0 = line.rstrip('\n')
    elif i == 1:
        t1 = line.rstrip('\n')
    elif i == 2:
        t2 = line.rstrip('\n')
    elif i == 3:
        t3 = line.rstrip('\n')
    i+=1

t0 = t0.split(',')
t1 = t1.split(',')
t2 = t2.split(',')
t3 = t3.split(',')

#x1= int(float(t0[0]))
#x2 =  int(float(t0[1]))
#x3 =  int(float(t0[2]))

#x4 = int(float(t0))
#date_FDM0 = datetime.date(x1 , x2,x3)
date_FDM0 = datetime.date(int(float(t0[0])),int(float(t0[1])),int(float(t0[2])))
date_FDMf = datetime.date(int(float(t1[0])),int(float(t1[1])),int(float(t1[2])))
date_IOP0 = datetime.date(int(float(t2[0])),int(float(t2[1])),int(float(t2[2])))
date_IOPf = datetime.date(int(float(t3[0])),int(float(t3[1])),int(float(t3[2])))


print 'FDM starts ' , date_FDM0 , ' and finshes ' , date_FDMf
print 'IOP starts ' , date_IOP0 , ' and finishes ', date_IOPf


#date_FDM0 = datetime.date(2017,3,21)
#date_FDMf = datetime.date(2017,3,28)
#date_IOP0 = datetime.date(2017,3,18)
#date_IOPf = datetime.date(2017,3,25)

#date_FDM0 = datetime.date(t0)
#date_FDMf = datetime.date(t1)
#date_IOP0 = datetime.date(t2)
#date_IOPf = datetime.date(t3)

n_days_FDM = (date_FDMf - date_FDM0).days + 1
n_days_IOP = (date_IOPf - date_IOP0).days + 1
dates_FDM = [(date_FDM0 + datetime.timedelta(n)).strftime('%Y%m%d') for n in range(n_days_FDM)]
dates_IOP = [(date_IOP0 + datetime.timedelta(n)).strftime('%Y%m%d') for n in range(n_days_IOP)]

def main():
    for date_FDM in dates_FDM:
        print "generating FDM report for " + date_FDM
        # process data using matlab function
        cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; daily_report_v6_nodisplay('%s','%s'); quit" """\
        % (matlab_path,myfun_path,"SIR_FDM_L2",date_FDM)
        cmd_FDM = shlex.split(cmd_FDM)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_FDM,stdout = devnull)
        
        # convert figures from eps to jpg
        cmd_convFig = "mogrify -density 300 -quality 100% -colorspace RGB -format jpg \
        /noc/users/cryo/QCV_Cryo2/code/gen_daily_report/figures/*FDM*.eps"
        cmd_convFig = shlex.split(cmd_convFig)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_convFig,stdout = devnull)
        
        # generate report
        genReportFun("FDM")

        # upload report to ESA server
        ftp_upDayReport("FDM",date_FDM)

    for date_IOP in dates_IOP:
        print "generating IOP report for " + date_IOP
        # process data using matlab function
        cmd_IOP = """%s -nodisplay -nosplash -nodesktop -r "cd %s; daily_report_v6_nodisplay('%s','%s'); quit" """\
        % (matlab_path,myfun_path,"SIR_IOP_L2",date_IOP)
        cmd_IOP = shlex.split(cmd_IOP)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_IOP,stdout = devnull)
        
        # convert figures from eps to jpg
        cmd_convFig = "mogrify -density 300 -quality 100% -colorspace RGB -format jpg \
        /noc/users/cryo/QCV_Cryo2/code/gen_daily_report/figures/*IOP*.eps"
        cmd_convFig = shlex.split(cmd_convFig)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_convFig,stdout = devnull)
        
        # generate report
        genReportFun("IOP")

        # upload report to ESA server
        ftp_upDayReport("IOP",date_IOP)
    

if __name__ == "__main__":
    main()
