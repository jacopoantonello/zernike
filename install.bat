@echo off
Powershell.exe -executionpolicy bypass -NoExit -Command ". scripts\base.ps1; Activate-Anaconda; pip uninstall -y zernike; git config core.autocrlf true; git config core.fileMode false; rm dist\*.whl -ErrorAction SilentlyContinue; python setup.py bdist_wheel; git config --unset core.autocrlf; git config --unset core.fileMode; pip install (get-item .\dist\*.whl);"



