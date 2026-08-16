@echo off
cd /d %~dp0
python campaign.py points ..\enrich_points.json ucs,bts,tx20 > enrich_console.log 2>&1
echo TERMINE >> enrich_console.log
