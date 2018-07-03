#./mc1_branching 500000 branching_250 600_harmonic.restlist 600_lower.restlist 40 40 250 &> branching_250.out &
#./mc1_branching 500000 branching_650 600_harmonic.restlist 600_lower.restlist 40 40 650 &> branching_650.out &

#test for non-empty branches
./mc1_branching_noempty 500000 branching_650 600_harmonic.restlist 600_lower.restlist 140 80 650 &> branching_650.out &
