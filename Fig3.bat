REM Fig 3 HR paper
REM vary NiBeta and q0
REM Fixed p0 = 0.3, Nialpha = 4, Nc = 0
REM
setlocal enableDelayedExpansion
@echo on
echo "Starting N=64, p0=0.3; Nc=0; Nialpha=4; NiBeta, q0 vary" !TIME!
REM
REM General batch file for Figure 3
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
REM Example call to Fig3.bat
REM
REM  Fig3.bat ring64\ output\ring64FixAvaryBNc.csv 1 10000
REM  Note:  For the paper, ran for size 64 spaces lattice, panmictic, ring,
REM         sfN64p1, sfN64p2, star
REM
fixed64 0.3 0.0 4 0 0.0 %4 %1 %3 1 > %2
REM
FOR %%b IN (2 4 8 16 32) DO (
fixed64 0.3 0.0 4 %%b 0.0 %4 %1 %3 0 >> %2
)

FOR %%b IN (0 2 4 8 16 32) DO (
FOR %%q IN (0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0) DO (
fixed64 0.3 %%q 4 %%b 0.0 %4 %1 %3 0 >> %2
)
)
