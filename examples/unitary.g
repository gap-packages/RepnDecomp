#! @BeginChunk Example_UnitaryRepresentation
#! @BeginExample
G := SymmetricGroup(3);;
irreps := IrreducibleRepresentations(G);;
# It happens that we are given unitary irreps, so
# rho is also unitary (its blocks are unitary)
rho := DirectSumOfRepresentations([irreps[1], irreps[2]]);;
IsUnitaryRepresentation(rho);
#! true
# Arbitrary change of basis
A := [ [ -1, 1 ], [ -2, -1 ] ];;
tau := ComposeHomFunction(rho, x -> A^-1 * x * A);;
# Not unitary, but still isomorphic to rho
IsUnitaryRepresentation(tau);
#! false
AreRepsIsomorphic(rho, tau);
#! true
# Now we unitarise tau
tau_u := UnitaryRepresentation(tau);;
# We get a record with the unitarised rep:
AreRepsIsomorphic(tau, tau_u.unitary_rep);
#! true
AreRepsIsomorphic(rho, tau_u.unitary_rep);
#! true
# The basis change is also in the record:
ForAll(G, g -> tau_u.basis_change * Image(tau_u.unitary_rep, g) = Image(tau, g) * tau_u.basis_change);
#! true
#! @EndExample
#! @EndChunk

#! @BeginChunk Example_IsUnitaryRepresentation
#! @BeginExample
# TODO: this example
#! @EndExample
#! @EndChunk


#! @BeginChunk Example_LDLDecomposition
#! @BeginExample
A := [ [ 3, 2*E(3)+E(3)^2, -3 ], [ E(3)+2*E(3)^2, -3, 3 ], [ -3, 3, -6 ] ];;
# A is a conjugate symmetric matrix
A = ConjugateTranspose@RepnDecomp(A);
#! true
# Note that A is not symmetric - the LDL decomposition works for any
# conjugate symmetric matrix.
A = TransposedMat(A);
#! false
decomp := LDLDecomposition(A);;
# The LDL decomposition is such that A = LDL^*, D diagonal, and L lower triangular.
A = decomp.L * DiagonalMat(decomp.D) * ConjugateTranspose@RepnDecomp(decomp.L);
#! true
decomp.L[1][2] = 0 and decomp.L[1][3] = 0 and decomp.L[2][3] = 0;
#! true
#! @EndExample
#! @EndChunk
