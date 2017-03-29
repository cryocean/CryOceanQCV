import os, shutil
import re
import uuid
import subprocess, shlex
from datetime import datetime, date, timedelta
from appy.pod.renderer import Renderer

product = "FDM"
template_odt = 'CryOcean-QCV_DayRep_' + product
path = '/noc/users/fmc1q07/QCV_Cryo2/code/gen_daily_report/'
image_dir = path + 'figures/'
odt_file_ext = '.odt'
Report_file_ext = '.pdf'
fid = open(path + 'figures/' + product + '_reportData.txt','r')
# Set standard variables
Version = 1.0
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
theor_ocean = fid.readline().strip()
ii = fid.readline()
ii = ii.split()
total_nrecord = ii[0]
total_precord = ii[1]
ii = fid.readline()
ii = ii.split()
ocean_nrecord = ii[0]
ocean_precord = ii[1]
ssh_nrecord = fid.readline().strip()
ssh_precord = fid.readline().strip()
ssh_ptheor = fid.readline().strip()
ssh_noise = fid.readline().strip()
ssh_noise1Hz = fid.readline().strip()
ssh_noiseDiff = fid.readline().strip()
ssh_nrecordNOC = fid.readline().strip()
ssh_precordNOC = fid.readline().strip()
ssh_ptheorNOC = fid.readline().strip()
swh_nrecord = fid.readline().strip()
swh_precord = fid.readline().strip()
swh_ptheor = fid.readline().strip()
swh_noise = fid.readline().strip()
swh_noise1Hz = fid.readline().strip()
swh_noiseDiff = fid.readline().strip()
swh_nrecordNOC = fid.readline().strip()
swh_precordNOC = fid.readline().strip()
swh_ptheorNOC = fid.readline().strip()
s0_nrecord = fid.readline().strip()
s0_precord = fid.readline().strip()
s0_ptheor = fid.readline().strip() 
s0_noise = fid.readline().strip()
s0_noise1Hz = fid.readline().strip()
s0_noiseDiff = fid.readline().strip()   
s0_nrecordNOC = fid.readline().strip()
s0_precordNOC = fid.readline().strip()
s0_ptheorNOC = fid.readline().strip() 
misp_nrecord = fid.readline().strip()
misp_precord = fid.readline().strip()
misp_ptheor = fid.readline().strip()
fid.close()

#Read warnings if any
nl = sum(1 for line in open(path + 'figures/' + product + "_warnings.txt")) #number of lines in file
warn = [];
if product in "IOP":
    lat_uni = " days"
else:
    lat_uni = " hours"

warn_type = ("latency_fail","latency_median_high","latency_mean_high","orbit_dropout","ocean_dropout",\
"large_orbit_bias")
if nl == 0:
    warn.append("No warnings")
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
    	        warn.append(warn_type[2] + ": " + warns[index] + lat_uni)
    	        index += 1
    	        if index > nl-1:
    		    break
    if index < nl-1:
        if warns[index] in warn_type[3]:
    	    index += 1
    	    while warns[index] not in warn_type:
    	        s1 = warns[index].split()
    	        warn.append(warn_type[3] + ": " + s1[1] + "% of records over ocean for orbit " +  s1[0])
    	        index += 1
    	        if index > nl-1:
    		    break
    if index < nl-1:
        if warns[index] in warn_type[4]:
    	    index += 1
    	    while warns[index] not in warn_type:
    	        warn.append(warn_type[4] + ": " + warns[index] + "% of records over ocean for " + Report_Date)
    	        index += 1
    	        if index > nl-1:
    		    break
    if index < nl-1:
        if warns[index] in warn_type[5]:
    	    index += 1
    	    while warns[index] not in warn_type:
    	        warn.append(warn_type[5] + " suspected for orbit " + warns[index])
    	        index += 1
    	        if index > nl-1:
    		    break
# Generate the Report
Report_fn = template_odt + Report_Date.replace("/","")
# call libreoffice in server mode (needed to convert to pdf)
accept_socket = '"socket,host=localhost,port=2002;urp;"'
soffi = "/usr/lib64/libreoffice/program/soffice --invisible --headless --accept=%s" % accept_socket
soffi = shlex.split(soffi)
pro = subprocess.Popen(soffi)
# run renderer
renderer = Renderer(path+template_odt+odt_file_ext, globals(), path + Report_fn+Report_file_ext,\
pythonWithUnoPath='/usr/lib64/libreoffice/program/pyuno')
renderer.run()
pro.terminate() #close libreoffice
