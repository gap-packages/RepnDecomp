gap> G := SymmetricGroup(4);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumRepList([irreps[1], irreps[2], irreps[3]]);;
gap> diag_info := BlockDiagonalRepresentationFast(rho);;
gap> rho := diag_info.diagonal_rep;;
gap> cc := ConjugacyClasses(G);;
gap> cent_basis := List(diag_info.centralizer_basis, BlockDiagonalMatrix);; # given in block form, need to convert
gap> cent_basis := List(cent_basis, m -> 1/Sqrt(Trace(m * m)) * m);; # need to normalize (already ortho)
gap> ForAll(cc, cl -> Sum(cl, s -> Image(rho, s)) = ClassSumCentralizer(rho, cl, cent_basis));
true