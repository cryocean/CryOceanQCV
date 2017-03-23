#!/usr/bin/env bash

sleep 5s
ndays=5
for ((i=1;i<=ndays;i++)) 
do  
    echo $i
    python "/noc/users/fmc1q07/Cryosat/get_data/ftp_syncCryoFast.py"
    #/usr/local/MATLAB/R2010b/bin/matlab -nojvm -nodisplay -nosplash -r \
    #"addpath('/home/kiko/Python_test/');sum_t; quit;"
    #python "/home/kiko/Python_test/ftp_cryosat3.py"
    sleep 10s
done
