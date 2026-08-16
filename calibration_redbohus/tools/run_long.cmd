@echo off
cd /d %~dp0
python campaign.py points ..\cal_long_points.json tx20,tx50,tx75,tx100 > long_console.log 2>&1
echo TERMINE >> long_console.log
