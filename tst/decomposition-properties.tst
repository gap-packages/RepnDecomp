# Property tests for exact decomposition and simultaneous block diagonalization.

gap> START_TEST("decomposition-properties.tst");

gap> G := SymmetricGroup(3);; irreps := IrreducibleRepresentations(G);;
gap> blockRep := DirectSumOfRepresentations([irreps[1], irreps[3], irreps[3]]);;
gap> A := REPN_TEST_UpperUnitriangularMat(5);;
gap> rho := REPN_TEST_ConjugateRepresentation(blockRep, A);;
gap> info := REPN_ComputeUsingSerre(rho : irreps := irreps);;
gap> REPN_TEST_IsBlockDiagonalization(rho, info);
true
gap> REPN_TEST_IsCollectedDecomposition(rho, Filtered(info.decomposition, x -> Length(x) > 0));
true
gap> List(Filtered(info.decomposition, x -> Length(x) > 0), Length) = [1,2];
true
gap> Sum(Flat(info.decomposition), summand -> Length(summand.basis)) = 5;
true
gap> Length(info.centralizer_basis) = 1^2 + 2^2;
true

gap> collected := IrreducibleDecompositionCollected(rho);;
gap> REPN_TEST_IsCollectedDecomposition(rho, collected);
true
gap> REPN_TEST_IsInvariantDecomposition(rho, IrreducibleDecomposition(rho));
true

gap> rhoK := REPN_TEST_ConjugateRepresentation(blockRep, A);;
gap> infoK := REPN_ComputeUsingSerre(rhoK : irreps := irreps, use_kronecker);;
gap> REPN_TEST_IsBlockDiagonalization(rhoK, infoK);
true
gap> REPN_TEST_IsCollectedDecomposition(rhoK, Filtered(infoK.decomposition, x -> Length(x) > 0));
true

gap> rhoAlt := REPN_TEST_ConjugateRepresentation(blockRep, A);;
gap> infoAlt := REPN_ComputeUsingMyMethod(rhoAlt : irreps := irreps);;
gap> REPN_TEST_IsBlockDiagonalization(rhoAlt, infoAlt);
true
gap> REPN_TEST_IsCollectedDecomposition(rhoAlt, infoAlt.decomposition);
true

gap> canonical := CanonicalDecomposition(REPN_TEST_ConjugateRepresentation(blockRep, A) : irreps := irreps);;
gap> SortedList(List(canonical, Dimension)) = [1,4];
true
gap> ForAll(canonical, V -> IsGInvariant@RepnDecomp(rho, V));
true
gap> RankMat(Concatenation(List(canonical, V -> List(Basis(V))))) = 5;
true

gap> T := Group(());; trivial3 := FuncToHom@RepnDecomp(T, g -> IdentityMat(3));;
gap> trivialDecomp := IrreducibleDecomposition(trivial3);;
gap> Length(trivialDecomp) = 3 and REPN_TEST_IsInvariantDecomposition(trivial3, trivialDecomp);
true

gap> H := SymmetricGroup(4);; hirreps := IrreducibleRepresentations(H);;
gap> multiplicityFree := DirectSumOfRepresentations([hirreps[1], hirreps[3], hirreps[5]]);;
gap> mfDecomp := IrreducibleDecomposition(multiplicityFree);;
gap> Length(mfDecomp) = 3 and REPN_TEST_IsInvariantDecomposition(multiplicityFree, mfDecomp);
true

gap> perm := ActionHomomorphism(G, [1..3]);;
gap> permCanonical := CanonicalDecomposition(perm : irreps := irreps);;
gap> SortedList(List(permCanonical, Dimension)) = [1,2];
true
gap> linearPerm := PermToLinearRep(perm);;
gap> ForAll(permCanonical, V -> IsGInvariant@RepnDecomp(linearPerm, V));
true

gap> STOP_TEST("decomposition-properties.tst", 1);
