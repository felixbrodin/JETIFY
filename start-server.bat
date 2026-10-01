@echo off
chcp 65001 >nul
title Jetify - local server
echo.
echo Starting local server for the Jetify app...
echo Then open http://localhost:8000 in your browser.
echo Press Ctrl+C in this window to stop the server.
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0start-server.ps1" 8000
