#! @BeginChunk Example_REPN_ComputeUsingMyMethod
#! @BeginExample
G := SymmetricGroup(4);;
irreps := IrreducibleRepresentations(G);;
rho := DirectSumOfRepresentations([irreps[4], irreps[5]]);;
# Jumble rho up a bit so it's not so easy for the library.
A := [ [ 3, -3, 2, -4, 0, 0 ], [ 4, 0, 1, -5, 1, 0 ], [ -3, -1, -2, 4, -1, -2 ],
       [ 4, -4, -1, 5, -3, -1 ], [ 3, -2, 1, 0, 0, 0 ], [ 4, 2, 4, -1, -2, 1 ] ];;
rho := ComposeHomFunction(rho, B -> A^-1 * B * A);;
# We've constructed rho from two degree 3 irreps, so there are a few
# things we can check for correctness:
decomposition := REPN_ComputeUsingMyMethod(rho);;
# Two distinct irreps, so the centralizer has dimension 2
Length(decomposition.centralizer_basis) = 2;
#! true
# Two distinct irreps i.e. two invariant subspaces
Length(decomposition.decomposition) = 2;
#! true
# All subspaces are dimension 3
ForAll(decomposition.decomposition, Vs -> Length(Vs) = 1 and Dimension(Vs[1]) = 3);
#! true
# And finally, check that the block diagonalized representation
# computed is actually isomorphic to rho:
AreRepsIsomorphic(rho, decomposition.diagonal_rep);
#! true
#! @EndExample
#! @EndChunk

#! @BeginChunk Example_REPN_ComputeUsingMyMethodCanonical
#! @BeginExample
# This is the same example as before, but splits into canonical
# summands internally. It gives exactly the same results, up to
# isomorphism.
other_decomposition := REPN_ComputeUsingMyMethodCanonical(rho);;
Length(other_decomposition.centralizer_basis) = 2;
#! true
Length(other_decomposition.decomposition) = 2;
#! true
ForAll(other_decomposition.decomposition, Vs -> Length(Vs) = 1 and Dimension(Vs[1]) = 3);
#! true
AreRepsIsomorphic(rho, other_decomposition.diagonal_rep);
#! true
#! @EndExample
#! @EndChunk
