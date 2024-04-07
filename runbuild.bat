@echo off
setlocal

call :runme call  app\omni.code.bat --ext-folder exts --enable coderage.io.spacemouse %* <NUL

goto :eof

:runme
%*
goto :eof
