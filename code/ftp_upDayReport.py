"""
ftp_upDayReport(data_type, date_report)

Uploads Cryosat2 daily reports for either FDM or IOP from the NOCs server
to the ESA server (Cryosat_GMQ@131.176.221.91).

INPUTS :    data_type : string denoting the Cryosat product. It can be
                        either "FDM" or "IOP". If data_type is not given
                        reports for both FDM and IOP are uploaded
            date_report : string denoting the date the report refer to in
                          the format yyyymmdd. If date_report is not given
                          date_report refers to present - 1 for FDM and to
                          present - 4 for IOP

OUTPUTS : This function does not have outputs  

author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
"""
import datetime
import ftplib
tout = 120
HOST = '131.176.221.91'
USER = 'Cryosat_GMQ'
PASSWD = 'cr10s4tGMQ'
rootR = '/MD0_DATA/Cryosat_GMQ/NOCS-FTP/QC-and-Validation-Reports/'
#rootL = '/scratch/general/cryosat/newReports/'

rootL = '/noc/mpoc/cryo/reports/daily/'


def ftp_upDayReport(data_type = None, date_report = None):
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
    if not data_type:
        data_type = ("FDM","IOP")
    else:
        if not isinstance(data_type,basestring):
            raise Exception("data_type must be a string")
        data_type = (data_type,)
    for d_type in data_type:
        if not date_report:
            if d_type == "FDM":
                date_report = datetime.date.today() - datetime.timedelta(days = 1)
            else:
                date_report = datetime.date.today() - datetime.timedelta(days = 4)
            date_report = date_report.strftime('%Y%m%d')
        try:
            ftp.cwd(rootR + d_type + "_daily") #set remote current work directory
        except:
            print "Could not set current work directory"
        listd = ftp.nlst()
        if date_report[0:4] not in listd: #if directory does not exist, make it (year)
            ftp.mkd(date_report[0:4])
        ftp.cwd(date_report[0:4])
        listd = ftp.nlst()
        if date_report[4:6] not in listd: #if directory does not exist, make it (month)
            ftp.mkd(date_report[4:6])
        ftp.cwd(date_report[4:6])
        
        filename = "CryOcean-QCV_DayRep_" + d_type + date_report + ".pdf"
        try:
            with open(rootL + filename,'r') as f:
                ftp.storbinary("STOR %s" %filename,f)
        except:
            print "File " + filename + " could not be uploaded"
        date_report = None
    ftp.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        ftp_upDayReport()
    elif len(sys.argv) == 2:
        ftp_upDayReport(sys.argv[1])
    elif len(sys.argv) == 3:
        ftp_upDayReport(sys.argv[1],sys.argv[2])
