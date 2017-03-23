import os, shutil
import re
import uuid
import subprocess, shlex
from datetime import datetime, date, timedelta
import time
from appy.pod.renderer import Renderer

def genReportFun():
    template_odt = 'CryOcean-QCV_MonthRep'
    #pathOut = '/scratch/general/cryosat/newReports_monthly/' changed 26/10/16
    pathOut = '/noc/mpoc/cryo/reports/monthly/'
    path = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/'
    image_dir = path + 'figures/'
    odt_file_ext = '.odt'
    Report_file_ext = '.pdf' #extension of the output file
    fid_FDM = open(path + 'figures/' + 'FDM' + '_reportData.txt','r')
    fid_IOP = open(path + 'figures/' + 'IOP' + '_reportData.txt','r')
    fid_GOP = open(path + 'figures/' + 'GOP' + '_reportData.txt','r')
    fidAll = [fid_FDM,fid_IOP,fid_GOP]
    prod = ['FDM','IOP','GOP']
    fld = {} #dictionary where data will be saved
    product = 'IOP' # added just for the warning section. Removed once section changed
    
    
    # Set standard variables
    Version = '1.1r1'
    Author = 'C. Banks'
    Checked_by = 'P. Cipollini'
    date_issued = datetime.utcnow().strftime('%d/%m/%y')
    
    ######## start reading data from file #############################
    for n,fid in enumerate(fidAll):
        # read latency
        latency = fid.readline()
        latency = latency.split()
            
        fld['median_latency' + prod[n]] = latency[0]
        fld['latencyLow' + prod[n]] = latency[1]
        fld['latencyUp' + prod[n]] = latency[2]

        # read data
        xVar = []
        for x in range(3):
            xVar.append(fid.readline())
        fld['start_time' + prod[n]] = xVar[0].strip()
        fld['stop_time' + prod[n]] = xVar[1].strip()
        orbits = xVar[2].split()
        fld['orbit0' + prod[n]] = orbits[0]
        fld['orbitf' + prod[n]] = orbits[1]
        Report_Date = fld['start_time' + prod[n]][0:11]
        Report_Date = Report_Date.strip()
        Report_Date = Report_Date.replace(" ","/")
        Report_Date = time.strptime(Report_Date,'%d/%m/%Y')
        Report_Date_fn = time.strftime('%b%Y',Report_Date)
        Report_Date = time.strftime('%B %Y',Report_Date)
        fld['theor_total' + prod[n]] = fid.readline().strip()
        fld['theor_ocean' + prod[n]] = fid.readline().strip()
        ii = fid.readline()
        ii = ii.split()
        fld['total_nrecord' + prod[n]] = ii[0]
        fld['total_precord' + prod[n]] = ii[1]
        ii = fid.readline()
        ii = ii.split()
        fld['ocean_nrecord' + prod[n]] = ii[0]
        fld['ocean_precord' + prod[n]] = ii[1]
        if prod[n] in ['IOP','GOP']:
            fld['rmsEdge_SSH' + prod[n]] = fid.readline().strip()
            fld['rmsOut_SSH' + prod[n]] = fid.readline().strip()
            fld['rmsEdge_SWH' + prod[n]] = fid.readline().strip()
            fld['rmsOut_SWH' + prod[n]] = fid.readline().strip()
            fld['rmsEdge_wsp' + prod[n]] = fid.readline().strip()
            fld['rmsOut_wsp' + prod[n]] = fid.readline().strip()
        fld['flagSLA' + prod[n]] = fid.readline().strip()
        fld['flagSWH' + prod[n]] = fid.readline().strip()
        fld['flagSigma0' + prod[n]] = fid.readline().strip()
        fld['flagWsp' + prod[n]] = fid.readline().strip()
        fld['flagMisp' + prod[n]] = fid.readline().strip()
        fld['editSLA' + prod[n]] = fid.readline().strip()
        fld['editSLAstd' + prod[n]] = fid.readline().strip()
        fld['editIB' + prod[n]] = fid.readline().strip()
        fld['BiasedOrbit' + prod[n]] = fid.readline().strip()
        fld['editWet' + prod[n]] = fid.readline().strip()
        fld['editDry' + prod[n]] = fid.readline().strip()
        fld['editIono' + prod[n]] = fid.readline().strip()
        fld['editSsb' + prod[n]] = fid.readline().strip()
        fld['editS0SLA' + prod[n]] = fid.readline().strip()
        fld['editS0stdSLA' + prod[n]] = fid.readline().strip()
        fld['editS0' + prod[n]] = fid.readline().strip()
        fld['editS0std' + prod[n]] = fid.readline().strip()
        fld['editSWH' + prod[n]] = fid.readline().strip()
        fld['editSWHstd' + prod[n]] = fid.readline().strip()
        fld['editWsp' + prod[n]] = fid.readline().strip()
        fld['editAllSLA' + prod[n]] = fid.readline().strip()
        fld['editAllSWH' + prod[n]] = fid.readline().strip()
        fld['editAllS0' + prod[n]] = fid.readline().strip()

        fid.close() # close file
    # convert dictionary to variables
    for key,val in fld.items():
        exec(key + '=val')

    #Read warnings if any
    nl = sum(1 for line in open(path + 'figures/' + product + "_warnings.txt")) #number of lines in file
    warn = [];
    if product in "IOP":
        lat_uni = " days"
    else:
        lat_uni = " hours"

    warn_type = ("latency_fail","latency_mean_high","orbit_dropout","ocean_dropout",\
    "large_orbit_bias")
    if nl == 0:
        warn = []
    else:
        warns = open(path + 'figures/' + product + "_warnings.txt").read().splitlines()
        index = 0
        if warns[index] in warn_type[0]:
            index += 1
            while warns[index] not in warn_type:
    	        warn.append(warn_type[0] + ": " + warns[index] + "% of records delivered within 3" + lat_uni)
    	        index += 1
    	        if index > nl-1:
    	            break
        if index < nl-1:
            if warns[index] in warn_type[1]:
    	        index += 1
    	        while warns[index] not in warn_type:
    	            warn.append(warn_type[1] + ": " + warns[index] + lat_uni)
    	            index += 1
    	            if index > nl-1:
    		        break
        if index < nl-1:
            if warns[index] in warn_type[2]:
    	        index += 1
    	        while warns[index] not in warn_type:
    	            s1 = warns[index].split()
    	            warn.append(warn_type[2] + ": " + s1[1] + "% of records over ocean for orbit " +  s1[0])
    	            index += 1
    	            if index > nl-1:
    		        break
        if index < nl-1:
            if warns[index] in warn_type[3]:
    	        index += 1
    	        while warns[index] not in warn_type:
    	            warn.append(warn_type[3] + ": " + warns[index] + "% of records over ocean for " + Report_Date)
    	            index += 1
    	            if index > nl-1:
    		        break
        if index < nl-1:
            if warns[index] in warn_type[4]:
    	        index += 1
    	        while warns[index] not in warn_type:
    	            warn.append(warn_type[4] + " suspected for orbit " + warns[index])
    	            index += 1
    	            if index > nl-1:
    		        break

    # Generate the Report
    Report_fn = template_odt + "_" + Report_Date_fn

    # call libreoffice in server mode (needed to convert to pdf)
    accept_socket = '"socket,host=localhost,port=2002;urp;"'
    soffi = "/usr/lib64/libreoffice/program/soffice --invisible --headless --accept=%s" % accept_socket
    soffi = shlex.split(soffi)
    pro = subprocess.Popen(soffi)

    # run renderer
    renderer = Renderer(path+template_odt+odt_file_ext, locals(), path + Report_fn+Report_file_ext,\
    pythonWithUnoPath='/usr/lib64/libreoffice/program/pyuno',overwriteExisting = True)
    renderer.run()
    pro.terminate() #close libreoffice
    
    # save file
    if os.path.isfile(pathOut + Report_fn + Report_file_ext):
        os.remove(pathOut + Report_fn + Report_file_ext)
    shutil.move(path + Report_fn + Report_file_ext,pathOut)    

if __name__ == "__main__":
    #import sys
    #genReportFun(sys.argv[1])
    genReportFun()
