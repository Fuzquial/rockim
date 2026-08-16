@echo off
cd /d %~dp0
python campaign.py points ..\cal_final_points.json ucs,bts,tx20,tx50 > final_console.log 2>&1
echo TERMINE >> final_console.log
