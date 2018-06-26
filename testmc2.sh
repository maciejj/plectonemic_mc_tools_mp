rm -f *test_run*
cp m10_40_40_600_250_mc2test.coord m10_40_40_600_250_mc2test_run.coord
cp m10_40_40_600_250_mc2test.pdb m10_40_40_600_250_mc2test_run.pdb

./mc2 20000 m10_40_40_600_250_mc2test_run 600_harmonic.restlist 600_lower.restlist 40 40 >& m10_40_40_600_250_mc2test_run.out
