import os, subprocess, shlex
import sys
import datetime
path_report = "/noc/users/cryo/QCV_Cryo2/code/gen_daily_report/"
sys.path.append(path_report)
from genReportFun import genReportFun

#matlab_path = "/noc/packages/linux_emt64/matlab/v2013a/bin/matlab" changed 24 Jan 2017

matlab_path = "/nerc/packages/matlab/2013a/bin/matlab"


myfun_path = "/noc/users/cryo/QCV_Cryo2/code/"
date_FDM = datetime.date.today() - datetime.timedelta(days = 1)
date_IOP = datetime.date.today() - datetime.timedelta(days = 4)
date_FDM = date_FDM.strftime('%Y%m%d')
date_IOP = date_IOP.strftime('%Y%m%d')

def main():
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; daily_report_v6_nodisplay('%s','%s'); quit" """\
    % (matlab_path,myfun_path,"SIR_FDM_L2",date_FDM)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    cmd_IOP = """%s -nodisplay -nosplash -nodesktop -r "cd %s; daily_report_v6_nodisplay('%s','%s'); quit" """\
    % (matlab_path,myfun_path,"SIR_IOP_L2",date_IOP)
    cmd_IOP = shlex.split(cmd_IOP)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_IOP,stdout = devnull)
     
    cmd_convFig = "mogrify -density 300 -quality 100% -colorspace RGB -format jpg \
    /noc/users/cryo/QCV_Cryo2/code/gen_daily_report/figures/*.eps"
    cmd_convFig = shlex.split(cmd_convFig)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_convFig,stdout = devnull)

    genReportFun("FDM")
    genReportFun("IOP")

if __name__ == "__main__":
    main()
