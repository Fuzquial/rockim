@echo off
cd /d %~dp0
python campaign.py points ..\pareto_points.json ucs,bts,tx20,tx50 > points_console.log 2>&1
echo TERMINE >> points_console.log
