gap> A := RandomInvertibleMat(5);;
gap> A := A + TransposedMat(A);;
gap> r := LDLDecomposition(A);;
gap> A = (r.L * DiagonalMat(r.D) * ConjugateTranspose@RepnDecomp(r.L));
true