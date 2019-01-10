gap> Read("tst/utils.g");;
gap> G := SymmetricGroup(4);;
gap> irreps := IrreducibleRepresentations(G);;
gap> # degrees of reps are 3, 2, 3
gap> my_irreps := [irreps[2], irreps[3], irreps[4]];;
gap> rho := DirectSumRepList(my_irreps);;
gap> # if we decompose, it should break into 3 bits with correct dimensions
gap> decomp := BlockDiagonalRepresentationFast(rho);;
gap> spaces := Flat(decomp.decomposition);;
gap> SortedList(List(spaces, V -> Dimension(V))) = SortedList(List(my_irreps, rep -> Trace(Image(rep, One(G)))));
true
gap> List(spaces, space -> IsGInvariant(rho, space)); # check the spaces are actually subrepresentations
[true, true, true]