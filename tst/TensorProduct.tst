gap> G := SymmetricGroup(4);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := Random(irreps);; # doesn't change between runs
gap> tau := Random(irreps);;
gap> tensor := TensorProductOfRepresentations(rho, tau);;
gap> kronecker := g -> KroneckerProduct(Image(rho, g), Image(tau, g));;
gap> A := RandomMat(DegreeOfRepresentation(rho), DegreeOfRepresentation(tau));;
gap> ForAll(GeneratorsOfGroup(G),
> g -> WrapMatrix@RepnDecomp(kronecker(g) * Flat(A), DegreeOfRepresentation(tau)) = Image(tensor, g) * A);
true