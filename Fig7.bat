REM Fig 7 HR paper
REM Fixed p,q, Nialpha, Vary NiBeta, Nc
REM Used to examine haplotype behaviour
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
REM Example call to Fig7.bat
REM
REM  Fig7.bat ring64\ output\fig9.N64.ring.csv 1 10000
REM  Note:  For the paper, ran for size 64 spaces lattice, panmictic, ring,
REM         sfN64p1, sfN64p2, star
REM
fixed64 0.1 0.1 8 2 0.0 %4 %1 %3 1 > %2
REM
FOR %%b IN (4 7 8) DO (
fixed64 0.1 0.1 8 %%b 0.0 %4 %1 %3 0 >> %2
)

FOR %%b IN (2 4 7 8) DO (
FOR %%c IN (0.0 0.0625 0.25 0.5 1.0 4.0 8.0) DO (
fixed64 0.1 0.1 8 %%b %%c %4 %1 %3 0 >> %2
)
)
