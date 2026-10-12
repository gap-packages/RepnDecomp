# Unitarization, Hermitian linear algebra, and orbital centralizers.

gap> START_TEST("unitary-permutation-properties.tst");

gap> G := SymmetricGroup(3);; irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[1],irreps[2]]);;
gap> IsUnitaryRepresentation(rho);
true
gap> A := [[-1,1],[-2,-1]];;
gap> tau := REPN_TEST_ConjugateRepresentation(rho, A);;
gap> IsUnitaryRepresentation(tau);
false
gap> result := UnitaryRepresentation(tau);;
gap> IsUnitaryRepresentation(result.unitary_rep);
true
gap> AreRepsIsomorphic(tau, result.unitary_rep);
true
gap> ForAll(GeneratorsOfGroup(G), g -> result.basis_change*Image(result.unitary_rep,g) = Image(tau,g)*result.basis_change);
true

gap> hermitian := [[4,1+E(4)],[1-E(4),3]];;
gap> ldl := LDLDecomposition(hermitian);;
gap> hermitian = ldl.L*DiagonalMat(ldl.D)*ConjugateTranspose@RepnDecomp(ldl.L);
true
gap> ForAll([1..Length(ldl.L)], i -> ForAll([1..Length(ldl.L)], j -> j<=i or ldl.L[i][j]=0));
true
gap> LDLDecomposition(IdentityMat(4)) = rec(L:=IdentityMat(4),D:=[1,1,1,1]);
true

gap> matrices := [[[1,0],[0,0]], [[0,0],[0,1]]];;
gap> orth := OrthonormalBasis@RepnDecomp(matrices);;
gap> IsOrthonormalSet(orth, InnerProduct@RepnDecomp);
true
gap> VectorSpace(Cyclotomics,orth) = VectorSpace(Cyclotomics,matrices);
true

gap> P := SymmetricGroup(3);; perm := ActionHomomorphism(P,[1..3]);;
gap> orbitalBasis := RepresentationCentralizerPermRep@RepnDecomp(perm);;
gap> Length(orbitalBasis) = 2;
true
gap> linearPerm := PermToLinearRep(perm);;
gap> ForAll(orbitalBasis, M -> REPN_TEST_CommutesWithRepresentation(linearPerm,M));
true
gap> Sum(orbitalBasis) = List([1..3], i -> [1,1,1]);
true
gap> RankMat(List(orbitalBasis,Flat)) = 2;
true

gap> STOP_TEST("unitary-permutation-properties.tst", 1);
