@echo off
rem lanceur rockim-studio (depuis la racine du repo)
cd /d "%~dp0"
python -m rockim_studio %*
