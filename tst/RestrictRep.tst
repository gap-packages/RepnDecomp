gap> G := SmallGroup(77, 1);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumRepList([irreps[1], irreps[2], irreps[3]]);;
gap> V := Cyclotomics^3;;
gap> B := Basis(V);;
gap> RestrictRep@RepnDecomp(rho, B) = rho;
true
gap> V := VectorSpace(Cyclotomics, [[0,1,0]], Zero(V));;
gap> B := Basis(V);;
gap> RestrictRep@RepnDecomp(rho, B) = irreps[2];
true
gap> RestrictRep@RepnDecomp(rho, B) = irreps[1];
false