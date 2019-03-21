gap> G := SmallGroup(32, 5);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[10], irreps[19], irreps[20], irreps[20]]);;
gap> # note deg irreps[10] = 1, deg for 19 and 20 is 2
gap> cent_basis := CentralizerBlocksOfRepresentation(rho);;
gap> # so dimension of centralizer should be sum of square of multiplcities
gap> Length(cent_basis) = 1 + 1 + 2^2;
true
gap> # sizes of blocks should be degree * multiplicity
gap> List(cent_basis[1], Length) = [1*1, 1*2, 2*2];
true