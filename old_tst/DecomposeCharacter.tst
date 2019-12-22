gap> G := SmallGroup(12, 2);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[1], irreps[3], irreps[5], irreps[5]]);;
gap> v := IrrVectorOfRepresentation@RepnDecomp(rho, List(irreps, CharacterOfRepresentation@RepnDecomp));
[ 1, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0 ]