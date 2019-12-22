gap> G := SmallGroup(56, 3);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[1], irreps[2], irreps[5]]);;
gap> # rho is already block diagonal, so shouldn't change
gap> BlockDiagonalRepresentation(rho) = rho;
true
gap> # now scramble the basis a bit and check if we find the original
gap> A := [ [ 1, 0, -1, -1 ], [ -1, 1, -1, 1 ], [ -2, -1, -2, 0 ], [ -1, 2, -2, -1 ] ];;
gap> tau := ComposeHomFunction(rho, x -> A^-1 * x * A);;
gap> BlockDiagonalRepresentation(tau) = rho;
true