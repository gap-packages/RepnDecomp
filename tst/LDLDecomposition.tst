gap> A := [ [ 1, 0, -1, -1, -1 ], [ 1, -1, 1, -2, -1 ], [ -2, 0, -1, 2, -2 ], [ -1, 2, -3, -1, 3 ], [ 0, -2, 1, -4, 0 ] ];; # TODO: randomness here
gap> A := A + TransposedMat(A);;
gap> r := LDLDecomposition(A);;
gap> A = (r.L * DiagonalMat(r.D) * ConjugateTranspose@RepnDecomp(r.L));
true