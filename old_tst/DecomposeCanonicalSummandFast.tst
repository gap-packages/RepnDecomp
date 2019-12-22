gap> G := DirectProduct(SymmetricGroup(3), SymmetricGroup(3));;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[2], irreps[2], irreps[2]]);;
gap> # so the canonical summand is just the whole space
gap> summand := Cyclotomics^3;;
gap> # and each axis is an irreducible G-invariant space, we should get
gap> # something resembling that (all G-inv, right number of summands)
gap> decomp := DecomposeCanonicalSummandAlternate@RepnDecomp(rho, irreps[2], summand);;
gap> ForAll(decomp, r -> IsGInvariant@RepnDecomp(rho, VectorSpace(Cyclotomics, r.basis)));
true
gap> Length(decomp);
3
gap> # something more complicated
gap> rho := DirectSumOfRepresentations([irreps[1], irreps[2], irreps[2], irreps[3]]);;
gap> decomp := DecomposeCanonicalSummandAlternate@RepnDecomp(rho, irreps[2], VectorSpace(Cyclotomics, [[0,1,0,0],[0,0,1,0]]));;
gap> ForAll(decomp, r -> IsGInvariant@RepnDecomp(rho, VectorSpace(Cyclotomics, r.basis)));
true
gap> Length(decomp);
2
gap> decomp := DecomposeCanonicalSummandAlternate@RepnDecomp(rho, irreps[3], VectorSpace(Cyclotomics, [[0,0,0,1]]));;
gap> ForAll(decomp, r -> IsGInvariant@RepnDecomp(rho, VectorSpace(Cyclotomics, r.basis)));
true
gap> Length(decomp);
1