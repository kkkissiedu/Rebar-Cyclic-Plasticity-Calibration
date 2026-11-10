@echo off
rem Launches the calibration app. Run from any folder; uses the "rebar-calib" conda environment if conda is installed.
setlocal
cd /d "%~dp0"
where conda >nul 2>&1 && call conda activate rebar-calib
python -m calibration_app.main
if errorlevel 1 (
  echo.
  echo The app did not start. Run:  python scripts\check_setup.py
  pause
)
