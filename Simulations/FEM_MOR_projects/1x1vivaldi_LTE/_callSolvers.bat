@echo off
REM Frequency of analysis
set freq=7.5e9
REM Searches automatically the *.hfss file (must be unique!!!) in the directory to obtain the name 
REM of the project
for /F "tokens=*" %%* in ('dir /b *.hfss') do call :Sub %%~n*
goto :eof

:Sub
echo Name of Project : %*.hfss
echo Chosen Frequency : %freq% Hz
echo.

REM Copies the HFSS model files (mesh, materials, boundary conditions) into the current directory
cd %*.hfssresults
cd *.results
cd *.cmesh
copy current.* ..\..\..\*
cd ..\..\..\
echo.
echo.

echo Please check if model files have been copied then 
echo press a button to call the Mesh Reader ...
echo 
pause
echo.

REM calling the MeshReader to create the LTE model from the HFSS model files
call ..\_femSolvers\ANST_MeshReader.exe %*
echo.
echo.
echo Press a button to call the Wave Solver ...
echo 
pause
echo.

set path=%path%;..\..\_femSolvers\

cd lte_fileset
REM Modify the *.mpara file to set the polynomials order (the last number in the SOLID definition)
REM remember that HFSS polynomial order begins with 0 while LTE's one with 1 (3D tent basis)
call notepad %*.mpara
REM Modify the *.seinfo file to obtain infos on the phases of excitation of the ports - these
REM informations are dump only if the EM_WaveSolver is called with the +singleEnded option
call notepad %*.seinfo

REM Calls the EM_WaveSolver in order to build the FEM matrix A and the right-hand sides b related to 
REM  the feeding ports - the solutions-to-fields operator will be dumped using the options
REM +fieldValues +fieldFunctional - the +directMatlab option relies on the mex compiled ParDiSo
call EM_WaveSolver.exe %* %freq% +directMatlab +fieldValues +fieldFunctional
REM call EM_WaveSolver.exe %* %freq% +singleEnded +directMatlab +fieldValues +fieldFunctional
echo 
echo.
pause