@echo off
cd /d %~dp0
python campaign.py screen > screen_console.log 2>&1
echo TERMINE >> screen_console.log
