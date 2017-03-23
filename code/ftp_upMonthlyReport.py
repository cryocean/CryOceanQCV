"""
ftp_upMonthlyReport(date_report)

Uploads Cryosat2 Monthly reports from the NOCs server to the ESA server (Cryosat_GMQ@131.176.221.91).

INPUTS :            date_report : string denoting the date the report refer to in
                    the format mmmyyyy (ex. "May2015").

OUTPUTS : This function does not have outputs  

author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
"""
import datetime
import ftplib
tout = 120
HOST = '131.176.221.91'
USER = 'Cryosat_GMQ'
PASSWD = 'cr10s4tGMQ'
rootR = '/MD0_DATA/Cryosat_GMQ/NOCS-FTP/QC-and-Validation-Reports/GOM_bimonthly/'
# rootL = '/scratch/general/cryosat/newReports_monthly/'
rootL = '/noc/users/cryo/QCV_Cryo2/reports/monthly'
def ftp_upMonthlyReport(date_report):
    print "connecting to ftp server...\n"
    n_attempts = 0
    n_attempts_max = 10
    while True:
        try:
            ftp = ftplib.FTP(HOST,timeout = tout)
            connected = True
        except(ftplib.all_errors,socket.timeout):
            print "timed out while trying to connect to server\n"
            if n_attempts >= n_attempts_max:
                raise Exception("max retry attempts reached. Exiting...")
            print "retrying connection to server in 5 minutes\n"
            time.sleep(300)
            connected = False
        if connected == True:
            break

    print "loging into ftp server...\n"
    try:
        ftp.login(USER,PASSWD)
    except:
        print "login failed\n"
    
    try:
        ftp.cwd(rootR) #set remote current work directory
    except:
        print "Could not set current work directory"

    listd = ftp.nlst()
    if date_report[-4:] not in listd: #if directory does not exist, make it (year)
        ftp.mkd(date_report[-4:])
    ftp.cwd(date_report[-4:])
    
    filename = "CryOcean-QCV_MonthRep_" + date_report + ".pdf"
    '''''print filename
    '''''
    try:
        with open(rootL + filename,'r') as f:
            ftp.storbinary("STOR %s" %filename,f)
    except:
        print "File " + filename + " could not be uploaded"
    ftp.close()

if __name__ == "__main__":
    import sys
    ftp_upMonthlyReport(sys.argv[1])
