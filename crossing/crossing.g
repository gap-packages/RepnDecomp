# Some algorithms to do with permutations that are needed to compute
# bounds on the crossing number of K_{m, n}

LoadPackage("RepnDecomp");

Read("utils.g");
Read("cr.g");

Drop := Drop@RepnDecomp;
Take := Take@RepnDecomp;

# Q(p, q) is the number of adjacent transpositions needed to get from
# p to q^-1. By this I don't mean by conjugation, I mean viewing p and
# q^-1 as cyclic lists and swapping adjacent elements without caring
# what they are.

# This is the full matrix of Q for S_n. Q is (n-1)! x (n-1)!, very
# large.  This probably won't be needed most of the time.
Qmatrix := function(n)
    local G, ncycles, big, i, j, gc, g, h, T, Q;

    # The code here is taken from the function agens in cr.g
    gc:=grp(n);
    g:=gc.g;
    ncycles := gc.c;
    h:=Subgroup(g,[g.1,g.2]);
    T:=TwoOrbitNumbers(g,[h]);
    StructureConstantsTwoOrbitsAlgebra(T);
    G:=NullGraph(g);
    AddEdgeOrbit(G,[1,1^g.3],h);

    # Takes indices into ncycles and tells you codistance between
    # ncycles[p] and ncycles[q]
    Q := function(p, q)
        # g.4 inverts cycles
        return Distance(G, p, q^g.4);
    end;

    # This is the full, big matrix of Q (not filled out yet)
    # WARNING: this could be extremely big
    big := NullMat(Factorial(n-1), Factorial(n-1));

    # Fill out the matrix really naively and stupidly
    for i in [1..Factorial(n-1)] do
        for j in [1..Factorial(n-1)] do
            big[i][j] := Q(i, j);
        od;
    od;

    # gives the matrix and the indexing list, so the matrix_ij is the
    # distance between index[i] and index[j]^-1
    return rec(matrix := big,
               index := ncycles);
end;

# gives the action of G as permutations of the list of m cycles
ActionPermRep := function(m)
    local G, mcycle, mcycles, action, ret;

    # This is the group whose action doesn't change Q(p, q)
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

    return ActionHomomorphism(G, mcycles, action);
end;

# gives the regular representation of the group action
ActionRegularRep := function(m)
    # We want the action to be represented as permutation matrices
    # Conjugating by any of these matrices fixes Q. This is the
    # representation we are block diagonalizing.
    return ConvertRhoIfNeeded@RepnDecomp(ActionPermRep(m));
end;

# The matrix Q is in the centralizer ring of action_hom - it commutes
# with all of the image matrices. At this point, we need Sage to
# compute the irreps of S_m, since Sage can give us integer matrices -
# much nicer than GAP's Cyclotomics matrices.

# Uses the irreps of S_m x S_2 to calculate the parameters for the SDP
# we need to solve
CalculateSDP := function(m, irreps)
    local block_diag_info, nice_basis, centralizer_basis, norm_cent_basis, d, centralizer, mult_param, param_matrices, i, j, k, nice_cent_basis, action_hom, cent_basis, action_perm, pairs, product, failmat, c, nice_change, nice_change_inv;

    Print("Computing group action: ");
    action_perm := ActionPermRep(m);
    Print("done\n");

    Print("Computing centralizer basis: ");
    # these are written in the standard basis
    centralizer_basis := RepresentationCentralizerPermRep@RepnDecomp(action_perm);
    Print("done\n");

    action_hom := PermToLinearRep(action_perm);

    # See https://homepages.cwi.nl/~lex/files/symm.pdf for the method we
    # now apply to get a smaller semidefinite program.

    # First, we block diagonalize action_hom (TODO: give cent_basis here)
    #block_diag_info := BlockDiagonalRepresentationFast(action_hom, irreps);
    #nice_basis := block_diag_info.basis;
    Print("Decomposing group action: ");
    nice_basis := BlockDiagonalBasis(action_hom);
    Print("done\n");

    nice_change := TransposedMat(nice_basis);
    nice_change_inv := nice_change^-1;

    # This is the basis for the centralizer, written in the nice
    # basis, called E_i in the paper. I convert to full matrices
    # here. (TODO: use sparse matrices here or something?). These
    # matrices must be orthogonal.
    #centralizer_basis := List(block_diag_info.centralizer_basis, BlockDiagonalMatrix);

    # for each E_i, there is some E_{i^*} = E_i^T. This is also the
    # case for our matrices. We record these pairs below. Possibly i =
    # i^* when you have a block on the diagonal
    Print("Calculating {i, i*} pairs: ");
    pairs := List([1..Length(centralizer_basis)], x -> fail);
    for i in [1..Length(pairs)] do
        for j in [1..Length(pairs)] do
            if centralizer_basis[i] = TransposedMat(centralizer_basis[j]) then
                # pairs[i] = i^*
                pairs[i] := j;
            fi;
        od;
    od;
    Print("done\n");

    Print("Normalizing bases: ");
    # We normalize the basis for the centralizer, each matrix is still
    # written in the standard basis. These are the B_i.
    norm_cent_basis := List(centralizer_basis, E -> (Sqrt(Trace(E*TransposedMat(E)))^-1)*E);

    # d is the dimension of the centralizer ring (sum of squares of
    # multiplicities of each irrep)
    d := Length(norm_cent_basis);

    # The centralizer ring itself
    centralizer := VectorSpace(Cyclotomics, norm_cent_basis);
    norm_cent_basis := Basis(centralizer, norm_cent_basis);

    # This is the same but in the block diag basis. All calculations
    # should be faster when done with these since the are block diagonal
    nice_cent_basis := List(norm_cent_basis, B -> nice_change_inv * B * nice_change);

    # normalise
    #nice_cent_basis := List(nice_cent_basis, B -> (1/Sqrt(InnerProduct@RepnDecomp(B, B))) * B);

    nice_cent_basis := Basis(VectorSpace(Cyclotomics, nice_cent_basis), nice_cent_basis);
    Print("done\n");

    # TODO: SPEED THIS UP BY A LOT!
    # The multiplication params are defined by B_i B_j = \sum_k \lambda_{i,j}^k B_k
    # The convention I use is that mult_param[i][j][k] = \lambda_{i,j}^k
    Print("Calculating (L_k)_ij: ");

    mult_param := NullMat(d, d);
    for i in [1..d] do
        for j in [1..d] do
            mult_param[i][j] := Coefficients(nice_cent_basis,
                                             nice_cent_basis[i]*nice_cent_basis[j]);
        od;
    od;

    Print("done\n");

    # The list of matrices (L_k)_{i,j} = \lambda_{k,j}^i
    param_matrices := List([1..d], k -> NullMat(d, d));
    for k in [1..d] do
        for i in [1..d] do
            for j in [1..d] do
                param_matrices[k][i][j] := mult_param[k][j][i];
            od;
        od;
    od;


    # # makes a square matrix of fails
    # failmat := function(n)
    #     return Replicate@RepnDecomp(Replicate@RepnDecomp(fail, n), n);
    # end;

    # # The list of matrices (L_k)_{i,j} = \lambda_{k,j}^i. Since B_k
    # # B_j = \sum_i \lambda_{k,j}^i B_i, can take inner product with
    # # B_i to get coefficient. There are also some tricks to fill out
    # # other coefficients due to symmetries etc.
    # param_matrices := List([1..d], k -> failmat(d));
    # for k in [1..d] do
    #     for j in [1..d] do
    #         product := norm_cent_basis[k]*norm_cent_basis[j];
    #         for i in [1..d] do
    #             if param_matrices[k][i][j] = fail then
    #                 param_matrices[k][i][j] := InnerProduct@RepnDecomp(product, norm_cent_basis[i]);
    #             fi;

    #             c := param_matrices[k][i][j];

    #             # these are some other coefficients we now know (see
    #             # paper for why this is true)
    #             param_matrices[pairs[i]][pairs[j]][k] := c;
    #             param_matrices[j][pairs[k]][pairs[i]] := c;
    #             param_matrices[pairs[j]][pairs[i]][pairs[k]] := c;
    #             param_matrices[pairs[k]][j][i] := c;
    #             param_matrices[i][k][pairs[j]] := c;
    #         od;
    #     od;
    # od;

    return rec(centralizer_basis := norm_cent_basis, # the B_i
               nice_basis := nice_basis, # basis that nicely block diagonalises everything
               nice_change := nice_change, # basis change matrix
               nice_change_inv := nice_change_inv, # nice_change^-1
               param_matrices := param_matrices, # the L_k
               pairs := pairs, # pairs[i] = i^*
              );
end;

CheckMatInCentralizer := function(m, mat)
    local action_hom, G;
    action_hom := ActionRegularRep(m);
    G := Source(action_hom);
    return ForAll(GeneratorsOfGroup(G),
                  g -> Image(action_hom, g)*mat = mat*Image(action_hom, g));
end;

# Now use Sage to solve the SDP problem from the paper to get \alpha_m
