@echo off
REM %1: project name

REM check if user entered an argument
if [%1]==[] goto error1
if exist %1 goto error2

REM Make CMake file
mkdir %1\src
echo PROJECT(%1) > %1\src\CMakeLists.txt
echo FIND_PACKAGE(ITK) >> %1\src\CMakeLists.txt
echo IF(ITK_FOUND) >> %1\src\CMakeLists.txt
echo 	INCLUDE(${ITK_USE_FILE}) >> %1\src\CMakeLists.txt
echo ELSE(ITK_FOUND) >> %1\src\CMakeLists.txt
echo MESSAGE(FATAL_ERROR >> %1\src\CMakeLists.txt
echo "ITK not found. Please set ITK_DIR.") >> %1\src\CMakeLists.txt
echo ENDIF(ITK_FOUND) >> %1\src\CMakeLists.txt
echo. >> %1\src\CMakeLists.txt
echo ADD_EXECUTABLE(%1 %1.cxx )  >> %1\src\CMakeLists.txt
echo. >> %1\src\CMakeLists.txt
echo TARGET_LINK_LIBRARIES(%1 ITKCommon ITKIO) >> %1\src\CMakeLists.txt

REM Make soure code
echo #include "itkImage.h" >> %1\src\%1.cxx
echo #include ^<iostream^> >> %1\src\%1.cxx
echo. >> %1\src\%1.cxx
echo int main() >> %1\src\%1.cxx
echo { >> %1\src\%1.cxx
echo 	system("pause"); >> %1\src\%1.cxx
echo 	return 0; >> %1\src\%1.cxx
echo } >> %1\src\%1.cxx

cd %1\src
cmake-gui
exit

goto end

:error1
echo Error: Please enter project name as argument
echo.
goto end

:error2
echo Error: Project already exists. Enter a new project name.
echo.

:end