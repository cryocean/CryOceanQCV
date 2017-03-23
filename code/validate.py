import os, subprocess, shlex
import sys
import datetime
path2Code = "/noc/users/cryo/QCV_Cryo2/code/validation/code_valid/"
sys.path.append(path2Code)


def validation(date_report):
    # date_report = "May-2014"
    
    path_figures = "/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/"
    myfun_path = "/noc/users/cryo/QCV_Cryo2/code/validation/code_valid/"
    #matlab_path = "/Applications/MATLAB_R2012a.app/bin/matlab" #Kiko's mac
    #matlab_path = "/noc/packages/linux_emt64/matlab/v2013a/bin/matlab" #Unix machine
    matlab_path = "/nerc/packages/matlab/2013a/bin/matlab" # changed 24 Jan 2017
    
    print "starting absolute tide gauge validation"
    "absolute tide gauge validation"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; absolute_val_tg('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Absolute tide gauge validation completed"
    
    print "starting relative tide gauge validation"    
    "relative tide gauge validation"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; relative_val_tg('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Relative tide gauge validation completed"
    
    print "starting SWH validation against buoys"    
    "SWH validation against buoys"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; swh_wsp_val_v2('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "SWH validation against buoys completed"
       
    print "starting wind validation against buoys"     
    "wind validation against buoys"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; cwsp_val_v2('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Wind validation against buoys completed"
        
    print "starting SWH validation against WW3 model"
    "SWH validation against WW3 model"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; swh_hist_ww3('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "SWH validation against WW3 model completed"
        
    print "starting validation against HF radar data"
    "validation against HF radar data"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; HF_val('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Validation against HF radar data completed"
        
    print "starting validation against OSCAR - Atlantic"
    "validation against OSCAR - Atlantic"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; Current_OSCAR_val('%s','%s'); quit" """\
    % (matlab_path,myfun_path,date_report,"Atlantic")
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Validation against OSCAR - Atlantic completed"
        
    print "starting validation against OSCAR - Pacific"
    "validation against OSCAR - Pacific"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; Current_OSCAR_val('%s','%s'); quit" """\
    % (matlab_path,myfun_path,date_report,"Pacific")
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Validation against OSCAR - Pacific completed"
        
    print "mapping SLA for validation against steric heights"
    "map SLA for validation against steric heights"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; map_sla('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "mapping SLA for validation against steric heights completed"
    
    print "starting computation of steric heights"
    "Read Argo data, compute steric heights, and save them in Matlab structure"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; dyn_argo_profiles; quit" """\
    % (matlab_path,myfun_path)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Computation of steric heights completed"
        
    print "starting validation against steric heights"
    "validation against steric heights"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; read_dyn_profiles('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Validation against steric heights completed"
        
    print "starting validation against J2 and C2 from Rads"
    "validation against Jason-2 and CryoSat-2 from Rads"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; J2_C2_rads_Val('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Validation against J2 and C2 from Rads completed"
        
    print "starting validation of GMSL"
    "GMSL comparison"
    cmd_FDM = """%s -nodisplay -nosplash -nodesktop -r "cd %s; GMSL_computation('%s'); quit" """\
    % (matlab_path,myfun_path,date_report)
    cmd_FDM = shlex.split(cmd_FDM)
    with open(os.devnull, 'w') as devnull:
        dum = subprocess.call(cmd_FDM,stdout = devnull)
    print "Validation of GMSL completed"
   
    
"""
matlab_path = "/noc/packages/linux_emt64/matlab/v2013a/bin/matlab"
myfun_path = "/noc/users/cryo/QCV_Cryo2/code/"
date_FDM = datetime.date.today() - datetime.timedelta(days = 1)
date_IOP = datetime.date.today() - datetime.timedelta(days = 4)
date_FDM = date_FDM.strftime('%Y%m%d')
date_IOP = date_IOP.strftime('%Y%m%d')
"""

if __name__ == "__main__":
    import sys
    validation(sys.argv[1])
