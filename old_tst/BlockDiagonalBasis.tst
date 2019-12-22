gap> G := SmallGroup(78, 5);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[2], irreps[3]]);;
gap> # it's already block diagonal so should be standard
gap> BlockDiagonalBasisOfRepresentation(rho);
[ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ]
gap> # it's hard to test this because there are many choices for this usually