import os, shutil
import re
import uuid
import subprocess, shlex
from datetime import datetime, date, timedelta
from appy.pod.renderer import Renderer

def genReportFunTDS(product):
    template_odt = 'CryOcean-QCV_DayRep_TDS_' + product
   # pathOut = '/scratch/general/cryosat/newReports/'
    pathOut = '/noc/mpoc/cryo/reports/daily/TDS'
    path = '/noc/users/cryo/QCV_Cryo2/code/gen_daily_report/'
    image_dir = path + 'figures/I_O_P_TDS'
    odt_file_ext = '.odt'
    Report_file_ext = '.pdf'
    fid = open(path + 'figures/I_O_P_TDS' + product + '_reportData.txt','r')
    # Set standard variables
    Version = '0.0r1'
    Author = 'C. Banks'
    Checked_by = 'P. Cipollini'
    date_issued = datetime.utcnow().strftime('%d/%m/%y')
    # read latency
    latency = fid.readline()
    latency = latency.split()
    median_latency = latency[0]
    latencyLow = latency[1]
    latencyUp = latency[2]

    # read data
    xVar = []
    for x in range(3):
        xVar.append(fid.readline())
    start_time = xVar[0].strip()
    stop_time = xVar[1].strip()
    orbits = xVar[2].split()
    orbit0 = orbits[0]
    orbitf = orbits[1]
    Report_Date = start_time[0:11]
    Report_Date = Report_Date.strip()
    Report_Date = Report_Date.replace(" ","/")
    theor_total = fid.readline().strip()
    theor_ocean = fid.readline().strip()
    ii = fid.readline()
    ii = ii.split()
    total_nrecord = ii[0]
    total_precord = ii[1]
    ii = fid.readline()
    ii = ii.split()
    ocean_nrecord = ii[0]
    ocean_precord = ii[1]
    if product is 'FDM':
        ssh_nrecord = fid.readline().strip()
        ssh_precord = fid.readline().strip()
        ssh_ptheor = fid.readline().strip()
        ssh_noise = fid.readline().strip()
        ssh_noise1Hz = fid.readline().strip()
        ssh_noiseDiff = fid.readline().strip()
        ssh_noiseNOC = fid.readline().strip()
        ssh_noiseNOC1Hz = fid.readline().strip()
        ssh_noiseDiffNOC = fid.readline().strip()
        ssh_nrecordNOC = fid.readline().strip()
        ssh_precordNOC = fid.readline().strip()
        ssh_ptheorNOC = fid.readline().strip()
        swh_nrecord = fid.readline().strip()
        swh_precord = fid.readline().strip()
        swh_ptheor = fid.readline().strip()
        swh_noise = fid.readline().strip()
        swh_noise1Hz = fid.readline().strip()
        swh_noiseDiff = fid.readline().strip()
        swh_noiseNOC = fid.readline().strip()
        swh_noiseNOC1Hz = fid.readline().strip()
        swh_noiseDiffNOC = fid.readline().strip()
        swh_nrecordNOC = fid.readline().strip()
        swh_precordNOC = fid.readline().strip()
        swh_ptheorNOC = fid.readline().strip()
        s0_nrecord = fid.readline().strip()
        s0_precord = fid.readline().strip()
        s0_ptheor = fid.readline().strip() 
        s0_noise = fid.readline().strip()
        s0_noise1Hz = fid.readline().strip()
        s0_noiseDiff = fid.readline().strip()
        s0_noiseNOC = fid.readline().strip()
        s0_noiseNOC1Hz = fid.readline().strip()
        s0_noiseDiffNOC = fid.readline().strip()   
        s0_nrecordNOC = fid.readline().strip()
        s0_precordNOC = fid.readline().strip()
        s0_ptheorNOC = fid.readline().strip()
        wind_nrecord = fid.readline().strip()
        wind_precord = fid.readline().strip()
        wind_ptheor = fid.readline().strip()
        wind_nrecordNOC = fid.readline().strip()
        wind_precordNOC = fid.readline().strip()
        wind_ptheorNOC = fid.readline().strip() 
        misp_nrecord = fid.readline().strip()
        misp_precord = fid.readline().strip()
        misp_ptheor = fid.readline().strip()
        misp_ptheorPos = fid.readline().strip()
    else:
        ssh_nrecord = fid.readline().strip()
        ssh_precord = fid.readline().strip()
        ssh_ptheor = fid.readline().strip()
        ssh_noise_LRM = fid.readline().strip()
        ssh_noise_PLRM = fid.readline().strip()
        ssh_noise1Hz_LRM = fid.readline().strip()
        ssh_noise1Hz_PLRM = fid.readline().strip()
        ssh_noiseDiff = fid.readline().strip()
        ssh_noiseNOC_LRM = fid.readline().strip()
        ssh_noiseNOC_PLRM = fid.readline().strip()
        ssh_noise1HzNOC_LRM = fid.readline().strip()
        ssh_noise1HzNOC_PLRM = fid.readline().strip()
        ssh_noiseDiffNOC = fid.readline().strip()
        ssh_nrecordNOC = fid.readline().strip()
        ssh_precordNOC = fid.readline().strip()
        ssh_ptheorNOC = fid.readline().strip()
        swh_nrecord = fid.readline().strip()
        swh_precord = fid.readline().strip()
        swh_ptheor = fid.readline().strip()
        swh_noise_LRM = fid.readline().strip()
        swh_noise_PLRM = fid.readline().strip()
        swh_noise1Hz_LRM = fid.readline().strip()
        swh_noise1Hz_PLRM = fid.readline().strip()
        swh_noiseDiff = fid.readline().strip()
        swh_noiseNOC_LRM = fid.readline().strip()
        swh_noiseNOC_PLRM = fid.readline().strip()
        swh_noise1HzNOC_LRM = fid.readline().strip()
        swh_noise1HzNOC_PLRM = fid.readline().strip()
        swh_noiseDiffNOC = fid.readline().strip()
        swh_nrecordNOC = fid.readline().strip()
        swh_precordNOC = fid.readline().strip()
        swh_ptheorNOC = fid.readline().strip()
        s0_nrecord = fid.readline().strip()
        s0_precord = fid.readline().strip()
        s0_ptheor = fid.readline().strip()
        s0_noise_LRM = fid.readline().strip()
        s0_noise_PLRM = fid.readline().strip()
        s0_noise1Hz_LRM = fid.readline().strip()
        s0_noise1Hz_PLRM = fid.readline().strip()
        s0_noiseDiff = fid.readline().strip()
        s0_noiseNOC_LRM = fid.readline().strip()
        s0_noiseNOC_PLRM = fid.readline().strip()
        s0_noise1HzNOC_LRM = fid.readline().strip()
        s0_noise1HzNOC_PLRM = fid.readline().strip()
        s0_noiseDiffNOC = fid.readline().strip()
        s0_nrecordNOC = fid.readline().strip()
        s0_precordNOC = fid.readline().strip()
        s0_ptheorNOC = fid.readline().strip()
        wind_nrecord = fid.readline().strip()
        wind_precord = fid.readline().strip()
        wind_ptheor = fid.readline().strip()
        wind_nrecordNOC = fid.readline().strip()
        wind_precordNOC = fid.readline().strip()
        wind_ptheorNOC = fid.readline().strip()
        misp_nrecord = fid.readline().strip()
        misp_precord = fid.readline().strip()
        misp_ptheor = fid.readline().strip()
        misp_ptheorPos = fid.readline().strip()
    editSLA = fid.readline().strip()
    editSLAstd = fid.readline().strip()
    editIB = fid.readline().strip()
    BiasedOrbit = fid.readline().strip()
    editWet = fid.readline().strip()
    editDry = fid.readline().strip()
    editIono = fid.readline().strip()
    editSsb = fid.readline().strip()
    editS0SLA = fid.readline().strip()
    editS0stdSLA = fid.readline().strip()
    editS0 = fid.readline().strip()
    editS0std = fid.readline().strip()
    editSWH = fid.readline().strip()
    editSWHstd = fid.readline().strip()
    editWsp = fid.readline().strip()
    editAllSLA = fid.readline().strip()
    editAllSWH = fid.readline().strip()
    editAllS0 = fid.readline().strip()

    fid.close()

    #Read warnings if any
    nl = sum(1 for line in open(path + 'figures/I_O_P_TDS' + product + "_warnings.txt")) #number of lines in file
    warn = [];
    if product in "IOP":
        lat_uni = " days"
        latMAD = "3"
    else:
        lat_uni = " hours"
        latMAD = "6.4"

    warn_type = ("latency_fail","latency_mean_high","orbit_dropout","ocean_dropout",\
    "large_orbit_bias","iono_missing","dry_missing","wet_missing","atm_missing","ssb_missing")
    if nl == 0:
        warn = []
    else:
        warns = open(path + 'figures/I_O_P_TDS' + product + "_warnings.txt").read().splitlines()
        index = 0
        if warns[index] in warn_type[0]:
            index += 1
            while warns[index] not in warn_type:
    	        warn.append(warn_type[0] + ": " + warns[index] + "% of records delivered with a delay of more than " + latMAD + lat_uni)
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

        if index < nl-1:
            if warns[index] in warn_type[5]:
                index += 1
                warn.append(warn_type[5] + ": ionospheric correction is missing")
        if index < nl-1:
            if warns[index] in warn_type[6]:
                index += 1
                warn.append(warn_type[6] + ": dry tropospheric correction is missing")
        if index < nl-1:
            if warns[index] in warn_type[7]:
                index += 1
                warn.append(warn_type[7] + ": wet tropospheric correction is missing")
        if index < nl-1:
            if warns[index] in warn_type[8]:
                index += 1
                warn.append(warn_type[8] + ": atmospheric correction is missing")
        if index < nl-1:
            if warns[index] in warn_type[9]:
                index += 1
                warn.append(warn_type[9] + ": sea state bias correction is missing")

    # Generate the Report
    Report_Date_fn = Report_Date.replace("/","")
    Report_Date_fn = Report_Date_fn[-4:] + Report_Date_fn[2:4] + Report_Date_fn[0:2]
    Report_fn = template_odt + Report_Date_fn
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
    if os.path.isfile(pathOut + Report_fn + Report_file_ext):
        os.remove(pathOut + Report_fn + Report_file_ext)
    shutil.move(path + Report_fn + Report_file_ext,pathOut)    

if __name__ == "__main__":
    import sys
    genReportFunTDSApril2017(sys.argv[1])
