REM Fig 5 HR paper
REM vary NiBeta and Nc
REM
setlocal enableDelayedExpansion
@echo on
echo "Starting N=64, Nialpha=8; p0=q0=0.1; NiBeta, Nc vary" !TIME!
REM
REM General batch file for Figure 5
REM Needs to be run for different networks
REM
REM Arguments to batch file:
REM 
REM %1 = networks path (folder path)
REM %2 = output file full path to save results 
REM %3 = Number of spaces (normally 1)
REM %4 = Number of independent runs (normally 10000)
REM 
REM Need to do first loop since we create the file initially
REM After this we can just do the multiple loops for each combination
REM
REM Example call to Fig5.bat
REM
REM  Fig5.bat ring64\ output\ring64FixAvaryBNc.csv 1 10000
REM  Note:  For the paper, ran for size 64 spaces lattice, panmictic, ring,
REM         sfN64p1, sfN64p2, star
REM
fixed64 0.1 0.1 8 1 0.0 %4 %1 %3 1 > %2
REM
FOR %%b IN (2 3 4 6 7 8 10 12 14 16 18 20 22 24 26 28 30 32 34) DO (
fixed64 0.1 0.1 8 %%b 0.0 %4 %1 %3 0 >> %2
)

FOR %%b IN (1 2 3 4 6 7 8 10 12 14 16 18 20 22 24 26 28 30 32 34) DO (
FOR %%c IN (0.25 1.0 2.0 4.0 8.0 16.0) DO (
fixed64 0.1 0.1 8 %%b %%c %4 %1 %3 0 >> %2
)
)
