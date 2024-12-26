@echo off
for /F "tokens=*" %%* in ('dir /b *.eps') do call :Sub %%~n*
goto :eof

:Sub
epstopdf %*.eps