#! @BeginChunk Example_CanonicalDecomposition
#! @BeginExample
# This is the trivial group
G := Group(());;
# The trivial group has only one representation per degree, so a
# degree d representation decomposes into a single canonical summand
# containing the whole space
rho := FuncToHom@RepnDecomp(G, g -> IdentityMat(3));;
canonical_summands_G := CanonicalDecomposition(rho);
#! [ ( Cyclotomics^3 ) ]
# More interesting example, S_3
H := SymmetricGroup(3);;
# The standard representation: a permutation to the corresponding
# permutation matrix.
tau := FuncToHom@RepnDecomp(H, h -> PermutationMat(h, 3));;
# Two canonical summands corresponding to the degree 2 and
# trivial irreps (in that order)
List(CanonicalDecomposition(tau), Dimension);
#! [ 2, 1 ]
#! @EndExample
#! @EndChunk

#! @BeginChunk Example_IrreducibleDecomposition
#! @BeginExample
# The trivial group has 1 irrep of degree 1, so rho decomposes into 3
# lines.
irred_decomp_G := IrreducibleDecomposition(rho);
#! [ rec( basis := [ [ 1, 0, 0 ] ] ), rec( basis := [ [ 0, 1, 0 ] ] ),
#!   rec( basis := [ [ 0, 0, 1 ] ] ) ]
# The spaces are returned in this format - explicitly keeping the
# basis - since this basis block diagonalises rho into the irreps,
# which are the smallest possible blocks. This is more obvious with
# H.
irred_decomp_H := IrreducibleDecomposition(tau);
#! [ rec( basis := [ [ 1, 1, 1 ] ] ),
#!   rec( basis := [ [ 1, E(3), E(3)^2 ], [ 1, E(3)^2, E(3) ] ] ) ]
# Using the basis vectors given there block diagonalises tau into
# the two blocks corresponding to the two irreps:
nice_basis := [ [ 1, 1, 1 ], [ 1, E(3), E(3)^2 ], [ 1, E(3)^2, E(3) ] ];;
tau_diag := ComposeHomFunction(tau, X -> nice_basis^-1 * X * nice_basis);
#! [ (1,2,3), (1,2) ] -> [ [ [ 1, 0, 0 ], [ 0, E(3), 0 ], [ 0, 0, E(3)^2 ] ],
#!   [ [ 1, 0, 0 ], [ 0, 0, E(3)^2 ], [ 0, E(3), 0 ] ] ]
#! @EndExample
#! @EndChunk

#! @BeginChunk Example_REPN_ComputeUsingSerre
#! @BeginExample
# Does the same thing we have done in the examples above, but all in
# one step, with as many subcomputations reused as possible
REPN_ComputeUsingSerre(tau);
#! rec( basis := [ [ 1, 1, 1 ], [ 1, E(3), E(3)^2 ], [ 1, E(3)^2, E(3) ] ],
#!   centralizer_basis := [ [ [ [ 1 ] ], [ [ 0, 0 ], [ 0, 0 ] ] ],
#!       [ [ [ 0 ] ], [ [ 1, 0 ], [ 0, 1 ] ] ] ],
#!   decomposition := [ [ rec( basis := [ [ 1, 1, 1 ] ] ) ], [  ],
#!       [ rec( basis := [ [ 1, E(3), E(3)^2 ], [ 1, E(3)^2, E(3) ] ] ) ] ],
#!   diagonal_rep := [ (1,2,3), (1,2) ] ->
#!     [ [ [ 1, 0, 0 ], [ 0, E(3), 0 ], [ 0, 0, E(3)^2 ] ],
#!       [ [ 1, 0, 0 ], [ 0, 0, E(3)^2 ], [ 0, E(3), 0 ] ] ] )
# You can also do the computation in parallel:
REPN_ComputeUsingSerre(tau : parallel);;
# Or specify the irreps if you have already computed them:
irreps_H := IrreducibleRepresentations(H);;
REPN_ComputeUsingSerre(tau : irreps := irreps_H);;
#! @EndExample
#! @EndChunk
