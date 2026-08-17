@echo off
cd /d %~dp0
set OMP_NUM_THREADS=14
.\rockim_chk.exe calibration_redbohus\configs\demo_cutting_graded.cfg calibration_redbohus\runs\demo_cut_graded > calibration_redbohus\tools\cutg_console.log 2>&1
echo TERMINE >> calibration_redbohus\tools\cutg_console.log
