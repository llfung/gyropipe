@echo off

echo Copying required files

set accessKey=\\icnas4.cc.ic.ac.uk\lsf212\ssh\ssh-vps.ppk

set CUR_YYYY=%date:~6,4%
set CUR_MM=%date:~3,2%
set CUR_DD=%date:~0,2%
set CUR_HH=%time:~0,2%
set CUR_NN=%time:~3,2%
set CUR_SS=%time:~6,2%
if "%CUR_HH:~0,1%" == " " set CUR_HH=0%CUR_HH:~1,1%
if "%CUR_NN:~0,1%" == " " set CUR_NN=0%CUR_NN:~1,1%
if "%CUR_SS:~0,1%" == " " set CUR_SS=0%CUR_SS:~1,1%
set uploadDirectory=gyro_%CUR_YYYY%%CUR_MM%%CUR_DD%_%CUR_HH%%CUR_NN%%CUR_SS%

echo %uploadDirectory%
echo "Uploading started"
REM This is actual all one command, however separated to separate lines for readability

winscp.com /command ^
 "option batch abort" ^
 "option reconnecttime off" ^
 "open scp://lsf212@login.hpc.ic.ac.uk -privatekey=%accessKey% -rawsettings SendBuf=0 Compression=1 -timeout=999" ^
 "cd ~"^
 "mkdir %uploadDirectory%"^
 "mkdir %uploadDirectory%/program"^
 "mkdir %uploadDirectory%/utils"^
 "cd %uploadDirectory%"^
 "put *.pbs ./"^
 "put state.cdf.in ./"^
 "put GTD_lib.cdf ./"^
 "put parallel.h ./"^
 "put Makefile ./"^
 "put num_core.sh ./"^
 "put .\program\* ./program/"^
 "put .\utils\* ./utils/"^
 "exit"

set cdstr=cd %uploadDirectory%
echo %cdstr% > shell_tmp.sh
type %CD%\shell_script.sh >> shell_tmp.sh

echo "Making files"

echo SSH to run
putty.exe -load "CX1" -m shell_tmp.sh

echo "RUN successful"

del shell_tmp.sh

REM exit
