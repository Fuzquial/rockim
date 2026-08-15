@echo off
cd /d %~dp0
python campaign.py lhs 44 > lhs_console.log 2>&1
echo TERMINE >> lhs_console.log
