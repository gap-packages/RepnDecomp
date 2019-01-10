gap> Read("tst/utils.g");;
gap> G := SymmetricGroup(4);;
gap> irreps := IrreducibleRepresentations(G);;
gap> my_irreps := [irreps[2], irreps[3], irreps[4]];;
gap> rho := DirectSumRepList(my_irreps);;
gap> # if we decompose, it should break into 3 bits with correct dimensions
gap> decomp := BlockDiagonalRepresentationFast(rho);;
gap> tau := decomp.diagonal_rep;;
gap> # check it gives an isomorphic rep
gap> AreRepsIsomorphic(rho, tau);
true
gap> spaces := Flat(decomp.decomposition);;
gap> # check dims of spaces and degrees of reps match up
gap> SortedList(List(spaces, Dimension)) = SortedList(List(my_irreps, DegreeOfRepresentation));
true
gap> List(spaces, space -> IsGInvariant(rho, space)); # check the spaces are actually subrepresentations
[true, true, true]