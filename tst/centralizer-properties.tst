# Centralizer and class-sum properties motivated by symmetry-reduced SDPs.

gap> START_TEST("centralizer-properties.tst");

gap> G := DihedralGroup(8);; irreps := IrreducibleRepresentations(G);;
gap> blockRep := DirectSumOfRepresentations([irreps[4], irreps[4], irreps[5]]);;
gap> A := REPN_TEST_UpperUnitriangularMat(4);;
gap> rho := REPN_TEST_ConjugateRepresentation(blockRep, A);;
gap> centralizer := CentralizerOfRepresentation(rho);;
gap> Length(centralizer) = 2^2 + 1^2;
true
gap> ForAll(centralizer, C -> REPN_TEST_CommutesWithRepresentation(rho, C));
true
gap> RankMat(List(centralizer, Flat)) = Length(centralizer);
true
gap> centSpace := VectorSpace(Cyclotomics, centralizer);;
gap> IdentityMat(4) in centSpace;
true
gap> ForAll(centralizer, X -> ForAll(centralizer, Y -> X*Y in centSpace));
true

gap> blockCentralizer := CentralizerBlocksOfRepresentation(rho);;
gap> Length(blockCentralizer) = 5;
true
gap> ForAll(blockCentralizer, blocks -> Length(blocks) = 2);
true
gap> List(blockCentralizer[1], Length) = [2,2];
true

gap> diagRep := BlockDiagonalRepresentation(rho);;
gap> expandedBlocks := List(blockCentralizer, BlockDiagonalMatrix@RepnDecomp);;
gap> ForAll(expandedBlocks, C -> REPN_TEST_CommutesWithRepresentation(diagRep, C));
true

gap> unitaryRep := blockRep;; unitaryCentralizer := CentralizerOfRepresentation(unitaryRep);;
gap> orthogonalCentralizer := OrthonormalBasis@RepnDecomp(unitaryCentralizer);;
gap> IsOrthonormalSet(orthogonalCentralizer, InnerProduct@RepnDecomp);
true
gap> ForAll(ConjugacyClasses(G), class -> ClassSumCentralizer(unitaryRep, class, orthogonalCentralizer) = Sum(class, g -> Image(unitaryRep,g)));
true
gap> ForAll(ConjugacyClasses(G), class -> ClassSumCentralizerNC(unitaryRep, class, orthogonalCentralizer) = Sum(class, g -> Image(unitaryRep,g)));
true

gap> sizes := [rec(dimension:=1,nblocks:=2), rec(dimension:=3,nblocks:=1)];;
gap> standardBlocks := SizesToBlocks(sizes);;
gap> Length(standardBlocks) = 5;
true
gap> ForAll(standardBlocks, blocks -> DimensionsMat(BlockDiagonalMatrix@RepnDecomp(blocks)) = [5,5]);
true
gap> RankMat(List(standardBlocks, blocks -> Flat(BlockDiagonalMatrix@RepnDecomp(blocks)))) = 5;
true

gap> STOP_TEST("centralizer-properties.tst", 1);
