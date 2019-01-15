LoadPackage("RepnDecomp");;

# a fairly large, non solvable group
G := SymmetricGroup(6);;

# we aren't benchmarking irrep list calculation, so precalculate it
irreps := IrreducibleRepresentations(G);;

# a "large" degree representation to decompose
#rho := DirectSumRepList([irreps[4], irreps[6], irreps[4], irreps[1], irreps[1]]);;
rho := irreps[4];
