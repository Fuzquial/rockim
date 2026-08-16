@echo off
cd /d %~dp0
python campaign.py points ..\rerun_points.json tx20,tx50 > rerun_console.log 2>&1
echo TERMINE >> rerun_console.log
