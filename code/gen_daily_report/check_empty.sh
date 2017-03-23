#!/usr/bin/env bash
ftp ftp://cryosat353:NoHtpQvL@science-pds.cryosat.esa.int
FILE=text.txt
if [[ -s $FILE ]] ; then
echo "$FILE has data."
else
echo "$FILE is empty."
fi ;

#! /bin/sh
HOST='servername'
USER='username'
PASSWD='password'
LOCAL_FILES='/local/dir'

wget ftp://$HOST/DIR01/* -nc --ftp-user=$USER --ftp-password=$PASSWD
wput --disable-tls --basename=$LOCAL_FILES/ $LOCAL_FILES/* ftp://$USER:$PASSWD@$HOST/DIR02/