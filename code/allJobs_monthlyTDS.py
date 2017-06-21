import os, subprocess, shlex
import sys, shutil
import datetime
path_report = ["/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/",\
"/noc/users/cryo/QCV_Cryo2/code/validation/","/noc/users/cryo/QCV_Cryo2/code/"]
sys.path.extend(path_report)
from genMonthlyReportFunTDS import genReportFun
from shutil import copyfile
from validate import validation
from ftp_upMonthlyReport import ftp_upMonthlyReport
matlab_path = "/nerc/packages/matlab/2013a/bin/matlab" # changed 24 Jan 2017
#matlab_path = "/noc/packages/linux_emt64/matlab/v2013a/bin/matlab"
myfun_path = "/noc/users/cryo/QCV_Cryo2/code/"
#configDateFile = "/noc/users/cryo/QCV_Cryo2/code/monthlyconfig.txt"  # defines start and stop dates


def main():
    t0 = 'Aug-2013'
    tf = 'Aug-2013'

    dates = [t0]
    ct = 0
    while dates[ct] != tf:
        dtmp = datetime.datetime.strptime(dates[ct],'%b-%Y') + datetime.timedelta(days=40)
        dates.append(dtmp.strftime('%b-%Y'))
        ct += 1
        
    for month in dates:    
        print "generating reports for " + month
        
        ''' 
        print "processing FDM data..."
        cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; monthly_report_IOP_FDM_v3_nodisplay('%s','%s'); quit" """\
        % (matlab_path,myfun_path,"SIR_FDM_L2",month)
        cmd_FDM = shlex.split(cmd_FDM)
        #print "processing FDM data Matlab PLUS finished"        
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_FDM,stdout = devnull)
        print "processing of FDM completed"
        
        print "processing IOP data..."
        cmd_IOP = """%s -nodisplay -nosplash -nodesktop -r "cd %s; monthly_report_IOP_FDM_v3_nodisplay('%s','%s'); quit" """\
        % (matlab_path,myfun_path,"SIR_IOP_L2",month)
        cmd_IOP = shlex.split(cmd_IOP)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_IOP,stdout = devnull)
        print "processing of IOP completed"
        '''
        
        print "processing GOP data..."
             
        cmd_GOP = """%s -nodisplay -nosplash -nodesktop -r "cd %s; monthly_report_IOP_FDM_v3_nodisplay_TDSv2('%s','%s'); quit" """\
        % (matlab_path,myfun_path,"SIR_GOP_L2",month)
        cmd_GOP = shlex.split(cmd_GOP)
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_GOP,stdout = devnull)
        
        print "processing of GOP completed"
        
        # REPLACE fig 6 with version run from Matlab direct as for some weird reason this does not work in nodisplay mode
        copyfile('/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/TDS/fig_6/GOP_Fig_6.eps',\
        '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/TDS/GOP_Fig_6.eps')

        print "converting figure from eps to png (could take a while)..."
        cmd_convFig = "mogrify -density 300 -quality 100% -colorspace RGB -format jpg \
        /noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/TDS/*.eps"
        cmd_convFig = shlex.split(cmd_convFig)
        #cmd_convFig = "mogrify -format png \
        with open(os.devnull, 'w') as devnull:
            dum = subprocess.call(cmd_convFig,stdout = devnull)
        print "conversion of figures completed"
        
        print "generating monthly report"
        genReportFun() # generate monthly report
        print "monthly report created" 
        
        
if __name__ == "__main__":
    main()
