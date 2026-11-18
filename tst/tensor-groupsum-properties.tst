# Tensor identities and exact BSGS group-sum comparisons.

gap> START_TEST("tensor-groupsum-properties.tst");

gap> G := SymmetricGroup(4);; irreps := IrreducibleRepresentations(G);;
gap> rho := irreps[4];; tau := irreps[3];;
gap> tensor := TensorProductOfRepresentations(rho, tau);;
gap> kron := KroneckerProductOfRepresentations(rho, tau);;
gap> testMatrix := [[1,2],[3,4],[5,6]];;
gap> ForAll(GeneratorsOfGroup(G), function(g) local explicit; explicit := KroneckerProduct(Image(rho,g),Image(tau,g)); return WrapMatrix@RepnDecomp(explicit*Flat(testMatrix),2) = Image(tensor,g)*testMatrix; end);
true
gap> ForAll(GeneratorsOfGroup(G), g -> Image(kron,g) = KroneckerProduct(Image(rho,g),Image(tau,g)));
true
gap> tensorChar := CharacterOfTensorProductOfRepresentations(tensor);;
gap> ForAll(GeneratorsOfGroup(G), g -> tensorChar(g) = Trace(Image(rho,g))*Trace(Image(tau,g)));
true

gap> gens := GeneratorsOfGroup(G);; x := Image(tensor,gens[1]);; y := Image(tensor,gens[2]);;
gap> (x*y)*testMatrix = x*(y*testMatrix);
true
gap> One(x)*testMatrix = testMatrix;
true
gap> x^-1*(x*testMatrix) = testMatrix;
true

gap> GroupSumBSGS(G, g -> Image(rho,g)) = Sum(G, g -> Image(rho,g));
true
gap> GroupSumBSGS(G, g -> Image(kron,g)) = Sum(G, g -> Image(kron,g));
true
gap> H := SmallGroup(12,3);; hrep := IrreducibleRepresentations(H)[4];;
gap> GroupSumBSGS(H, g -> Image(hrep,g)) = Sum(H, g -> Image(hrep,g));
true
gap> T := Group(());; trep := FuncToHom@RepnDecomp(T, g -> [[1]]);;
gap> GroupSumBSGS(T, g -> Image(trep,g)) = [[1]];
true

gap> direct := DirectSumOfRepresentations([irreps[1],irreps[3]]);;
gap> ForAll(GeneratorsOfGroup(G), g -> Image(direct,g) = BlockDiagonalMatrix@RepnDecomp([Image(irreps[1],g),Image(irreps[3],g)]));
true

gap> G2 := SymmetricGroup(2);; G3 := SymmetricGroup(3);;
gap> products := TensorProductRepLists(IrreducibleRepresentations(G2), IrreducibleRepresentations(G3));;
gap> Length(products) = Length(IrreducibleRepresentations(G2))*Length(IrreducibleRepresentations(G3));
true
gap> ForAll(products, rep -> Source(rep) = DirectProduct(G2,G3));
true

gap> STOP_TEST("tensor-groupsum-properties.tst", 1);
