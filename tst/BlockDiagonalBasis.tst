gap> Read("tst/utils.g");;
gap> G := SmallGroup(78, 5);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumRepList([irreps[2], irreps[3]]);;
gap> # it's already block diagonal so should be standard
gap> BlockDiagonalBasis(rho);
[ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ]
gap> # it's hard to test this because there are many choices for this usually