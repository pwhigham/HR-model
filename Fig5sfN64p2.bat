REM Fig 2 HR paper
REM
setlocal enableDelayedExpansion
@echo on
echo "Starting N=64, Nialpha=8; NiBeta vary" !TIME!
REM p,q,Nialpha,Nibeta, Nc runs
fixed64 0.100 0.1 8 1 0.0 10000 networks\sf64p2\ 100 1 > output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 2 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 3 " !TIME!
fixed64 0.100 0.1 8 3 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 4 " !TIME!
fixed64 0.100 0.1 8 4 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 6 " !TIME!
fixed64 0.100 0.1 8 6 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 7 " !TIME!
fixed64 0.100 0.1 8 7 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 8 " !TIME!
fixed64 0.100 0.1 8 8 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 10 " !TIME!
fixed64 0.100 0.1 8 10 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 12" !TIME!
fixed64 0.100 0.1 8 12 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 14 " !TIME!
fixed64 0.100 0.1 8 14 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 16 " !TIME!
fixed64 0.100 0.1 8 16 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 18 " !TIME!
fixed64 0.100 0.1 8 18 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 20 " !TIME!
fixed64 0.100 0.1 8 20 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 22 " !TIME!
fixed64 0.100 0.1 8 22 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 24 " !TIME!
fixed64 0.100 0.1 8 24 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 26 " !TIME!
fixed64 0.100 0.1 8 26 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 28 " !TIME!
fixed64 0.100 0.1 8 28 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 30 " !TIME!
fixed64 0.100 0.1 8 30 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 32 " !TIME!
fixed64 0.100 0.1 8 32 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 34 0.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
REM p,q,Nialpha,Nibeta, Nc=0.25 runs
fixed64 0.100 0.1 8 1 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 2 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 3 " !TIME!
fixed64 0.100 0.1 8 3 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 4 " !TIME!
fixed64 0.100 0.1 8 4 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 6 " !TIME!
fixed64 0.100 0.1 8 6 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 7 " !TIME!
fixed64 0.100 0.1 8 7 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 8 " !TIME!
fixed64 0.100 0.1 8 8 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 10 " !TIME!
fixed64 0.100 0.1 8 10 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 12" !TIME!
fixed64 0.100 0.1 8 12 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 14 " !TIME!
fixed64 0.100 0.1 8 14 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 16 " !TIME!
fixed64 0.100 0.1 8 16 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 18 " !TIME!
fixed64 0.100 0.1 8 18 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 20 " !TIME!
fixed64 0.100 0.1 8 20 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 22 " !TIME!
fixed64 0.100 0.1 8 22 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 24 " !TIME!
fixed64 0.100 0.1 8 24 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 26 " !TIME!
fixed64 0.100 0.1 8 26 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 28 " !TIME!
fixed64 0.100 0.1 8 28 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 30 " !TIME!
fixed64 0.100 0.1 8 30 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 32 " !TIME!
fixed64 0.100 0.1 8 32 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 34 0.25 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
REM p,q,Nialpha,Nibeta, Nc runs
fixed64 0.100 0.1 8 1 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 2 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 3 " !TIME!
fixed64 0.100 0.1 8 3 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 4 " !TIME!
fixed64 0.100 0.1 8 4 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 6 " !TIME!
fixed64 0.100 0.1 8 6 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 7 " !TIME!
fixed64 0.100 0.1 8 7 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 8 " !TIME!
fixed64 0.100 0.1 8 8 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 10 " !TIME!
fixed64 0.100 0.1 8 10 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 12" !TIME!
fixed64 0.100 0.1 8 12 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 14 " !TIME!
fixed64 0.100 0.1 8 14 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 16 " !TIME!
fixed64 0.100 0.1 8 16 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 18 " !TIME!
fixed64 0.100 0.1 8 18 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 20 " !TIME!
fixed64 0.100 0.1 8 20 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 22 " !TIME!
fixed64 0.100 0.1 8 22 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 24 " !TIME!
fixed64 0.100 0.1 8 24 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 26 " !TIME!
fixed64 0.100 0.1 8 26 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 28 " !TIME!
fixed64 0.100 0.1 8 28 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 30 " !TIME!
fixed64 0.100 0.1 8 30 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 32 " !TIME!
fixed64 0.100 0.1 8 32 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 34 1.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
REM p,q,Nialpha,Nibeta, Nc runs
fixed64 0.100 0.1 8 1 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 2 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 3 " !TIME!
fixed64 0.100 0.1 8 3 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 4 " !TIME!
fixed64 0.100 0.1 8 4 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 6 " !TIME!
fixed64 0.100 0.1 8 6 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 7 " !TIME!
fixed64 0.100 0.1 8 7 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 8 " !TIME!
fixed64 0.100 0.1 8 8 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 10 " !TIME!
fixed64 0.100 0.1 8 10 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 12" !TIME!
fixed64 0.100 0.1 8 12 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 14 " !TIME!
fixed64 0.100 0.1 8 14 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 16 " !TIME!
fixed64 0.100 0.1 8 16 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 18 " !TIME!
fixed64 0.100 0.1 8 18 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 20 " !TIME!
fixed64 0.100 0.1 8 20 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 22 " !TIME!
fixed64 0.100 0.1 8 22 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 24 " !TIME!
fixed64 0.100 0.1 8 24 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 26 " !TIME!
fixed64 0.100 0.1 8 26 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 28 " !TIME!
fixed64 0.100 0.1 8 28 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 30 " !TIME!
fixed64 0.100 0.1 8 30 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 32 " !TIME!
fixed64 0.100 0.1 8 32 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 34 4.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
REM p,q,Nialpha,Nibeta, Nc runs
fixed64 0.100 0.1 8 1 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 2 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 3 " !TIME!
fixed64 0.100 0.1 8 3 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 4 " !TIME!
fixed64 0.100 0.1 8 4 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 6 " !TIME!
fixed64 0.100 0.1 8 6 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 7 " !TIME!
fixed64 0.100 0.1 8 7 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 8 " !TIME!
fixed64 0.100 0.1 8 8 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 10 " !TIME!
fixed64 0.100 0.1 8 10 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 12" !TIME!
fixed64 0.100 0.1 8 12 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 14 " !TIME!
fixed64 0.100 0.1 8 14 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 16 " !TIME!
fixed64 0.100 0.1 8 16 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 18 " !TIME!
fixed64 0.100 0.1 8 18 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 20 " !TIME!
fixed64 0.100 0.1 8 20 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 22 " !TIME!
fixed64 0.100 0.1 8 22 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 24 " !TIME!
fixed64 0.100 0.1 8 24 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 26 " !TIME!
fixed64 0.100 0.1 8 26 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 28 " !TIME!
fixed64 0.100 0.1 8 28 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 30 " !TIME!
fixed64 0.100 0.1 8 30 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 32 " !TIME!
fixed64 0.100 0.1 8 32 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 34 8.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
REM p,q,Nialpha,Nibeta, Nc runs
fixed64 0.100 0.1 8 1 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 2 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 3 " !TIME!
fixed64 0.100 0.1 8 3 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 4 " !TIME!
fixed64 0.100 0.1 8 4 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 6 " !TIME!
fixed64 0.100 0.1 8 6 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 7 " !TIME!
fixed64 0.100 0.1 8 7 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 8 " !TIME!
fixed64 0.100 0.1 8 8 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 10 " !TIME!
fixed64 0.100 0.1 8 10 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 12" !TIME!
fixed64 0.100 0.1 8 12 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 14 " !TIME!
fixed64 0.100 0.1 8 14 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 16 " !TIME!
fixed64 0.100 0.1 8 16 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 18 " !TIME!
fixed64 0.100 0.1 8 18 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 20 " !TIME!
fixed64 0.100 0.1 8 20 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 22 " !TIME!
fixed64 0.100 0.1 8 22 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 24 " !TIME!
fixed64 0.100 0.1 8 24 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 26 " !TIME!
fixed64 0.100 0.1 8 26 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 28 " !TIME!
fixed64 0.100 0.1 8 28 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 30 " !TIME!
fixed64 0.100 0.1 8 30 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 32 " !TIME!
fixed64 0.100 0.1 8 32 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
fixed64 0.100 0.1 8 34 16.0 10000 networks\sf64p2\ 100 0 >> output\sf64p2FixAvaryBNc.csv
echo "Beta 2 " !TIME!
