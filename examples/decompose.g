LoadPackage("RepnDecomp");

G := SymmetricGroup(3);

# f is the permutation representation of the action of G on itself
f := GroupHomomorphismByImages(G, Range(RegularActionHomomorphism(G)));

# rho is the linear representation, which is what our package needs
rho := PermToLinearRep(f);

# We can now decompose it into irreducibles
Print(IrreducibleDecomposition(rho));
