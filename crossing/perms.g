# Some algorithms to do with permutations that are needed to compute
# bounds on the crossing number of K_{m, n}

LoadPackage("RepnDecomp");

Read("compute_q.g");

Drop := Drop@RepnDecomp;
Take := Take@RepnDecomp;

# Q(p, q) is the number of adjacent transpositions needed to get from
# p to q^-1. By this I don't mean by conjugation, I mean viewing p and
# q^-1 as cyclic lists and swapping adjacent elements without caring
# what they are.
#
# NOTE: This is not a clever way to calculate this, but it works fine.
Q := function(p, q) return AdjacentTranspositionsBetween(p, q^(-1)); end;

# This is the full matrix of Q for S_n. Q is (n-1)! x (n-1)!
# This probably won't be needed most of the time
Qmatrix := function(n)
    local G, ncycle, ncycles, big, i, j;

    G := SymmetricGroup(n);

    # This is just some n-cycle
    ncycle := MappingPermListList(ShiftLeft([1..n], 1), [1..n]);

    # The n-cycles index the rows and cols of Q
    ncycles := List(ncycle^G);

    # This is the full, big matrix of Q (not filled out yet)
    # WARNING: this could be extremely big
    big := NullMat(Factorial(n-1), Factorial(n-1));

    # Fill out the matrix really naively and stupidly
    for i in [1..Factorial(n-1)] do
        for j in [1..Factorial(n-1)] do
            big[i][j] := AdjacentTranspositionsBetween(ncycles[i], ncycles[j]^-1);
        od;
    od;

    # gives the matrix and the indexing list, so Qmatrix[i][j] =
    # Q(index[i], index[j]).
    return rec(matrix := big,
               index := ncycles);
end;

# Here we fix some number, the end goal is to calculate \alpha_m
m := 3;

# This is the group Q is invariant under
G := DirectProduct(SymmetricGroup(m), SymmetricGroup(2));

# a random m cycle
mcycle := MappingPermListList(ShiftLeft([1..m], 1), [1..m]);

# m-cycles index the rows and cols of Q
mcycles := List(mcycle^G);

# This is how G acts on mcycles. h_(pi, i)(p) = pi p^i pi^-1
action := function(cycle, g)
    local g1, g2, i, real_i, pi;

    # Acts via conjugation
    pi := Image(Projection(G, 1), g);

    # Acts by inverting the cycle
    i := Image(Projection(G, 2), g);

    real_i := 1;

    # if it's nontrivial, invert
    if i <> One(SymmetricGroup(2)) then
        real_i := -1;
    fi;

    return pi * cycle^real_i * (pi^-1);
end;

action_hom := ActionHomomorphism(G, mcycles, action);

# We want the action to be represented as permutation matrices
# Conjugating by any of these matrices fixes Q. This is the
# representation we are block diagonalizing.
action_hom := ConvertRhoIfNeeded@RepnDecomp(action_hom);

# The matrix Q is in the centralizer ring of action_hom - it commutes
# with all of the image matrices. At this point, we need Sage to
# compute the irreps of S_m, since Sage can give us integer matrices -
# much nicer than GAP's Cyclotomics matrices.

# Uses the irreps of S_m x S_2 to calculate the parameters for the SDP we need to solve
CalculateSDP := function(irreps)
    local block_diag_info, nice_basis, centralizer_basis, norm_cent_basis, d, centralizer, mult_param, param_matrices, i, j, k, nice_cent_basis;

    # See https://homepages.cwi.nl/~lex/files/symm.pdf for the method we
    # now apply to get a smaller semidefinite program.

    # First, we block diagonalize action_hom
    block_diag_info := BlockDiagonalRepresentationFast(action_hom, irreps);
    nice_basis := block_diag_info.basis;

    # This is the nice basis for the centralizer, written in the nice
    # basis, called E_i in the paper. I convert to full matrices
    # here. (TODO: use sparse matrices here or something?)
    centralizer_basis := List(block_diag_info.centralizer_basis, blocks -> BlockDiagonalMatrix(blocks));

    # We normalize the basis for the centralizer, each matrix is still
    # written in the nice basis. These are the B_i.
    norm_cent_basis := List(centralizer_basis, E -> (Sqrt(Trace(TransposedMat(E)*E))^-1)*E);

    # d is the dimension of the centralizer ring (sum of squares of
    # multiplicities of each irrep)
    d := Length(norm_cent_basis);

    # The centralizer ring itself
    centralizer := VectorSpace(Cyclotomics, norm_cent_basis);
    nice_cent_basis := Basis(centralizer, norm_cent_basis);

    # The multiplication params are defined by B_i B_j = \sum_k \lambda_{i,j}^k B_k
    # The convention I use is that mult_param[i][j][k] = \lambda_{i,j}^k
    mult_param := NullMat(d, d);
    for i in [1..d] do
        for j in [1..d] do
            mult_param[i][j] := Coefficients(nice_cent_basis,
                                             norm_cent_basis[i]*norm_cent_basis[j]);
        od;
    od;

    # The list of matrices (L_k)_{i,j} = \lambda_{k,j}^i
    param_matrices := List([1..d], k -> NullMat(d, d));
    for k in [1..d] do
        for i in [1..d] do
            for j in [1..d] do
                param_matrices[k][i][j] := mult_param[k][j][i];
            od;
        od;
    od;

    return rec(centralizer_basis := norm_cent_basis, # the B_i
               nice_basis := nice_basis, # basis all matrices are written in
               mult_param := mult_param, # the lambda_{i,j}^k
               param_matrices := param_matrices); # the L_k
end;

# Now use Sage to solve the SDP problem from the paper to get \alpha_m
