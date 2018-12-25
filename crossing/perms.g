# Some algorithms to do with permutations that are needed to compute
# bounds on the crossing number of K_{m, n}

LoadPackage("RepnDecomp");

Read("compute_q.g");

Drop := Drop@RepnDecomp;
Take := Take@RepnDecomp;

# swaps l[pos] left one place (cyclically)
swap_left := function(l, pos)
    local tmp;
    if pos = 1 then
        tmp := l[pos];
        l[pos] := l[Length(l)];
        l[Length(l)] := tmp;
    else
        tmp := l[pos];
        l[pos] := l[pos-1];
        l[pos-1] := tmp;
    fi;
end;

# swaps l[pos] right one place (cyclically)
swap_right := function(l, pos)
    local tmp;
    if pos = Length(l) then
        tmp := l[pos];
        l[pos] := l[1];
        l[1] := tmp;
    else
        tmp := l[pos];
        l[pos] := l[pos+1];
        l[pos+1] := tmp;
    fi;
end;

# Number of adjacent swaps needed to sort a permutation of [1..n],
# where we consider the first and last elements to be adjacent. We
# don't move the element "fixed".
SwapsToSortWithFixed := function(l, fixed)
    local list, n, nswaps, pos_fixed, elem, curpos, correct_offset, current_offset, left_list, left_moves, right_list, right_moves;
    list := ShallowCopy(l);
    n := Maximum(list);

    nswaps := 0;

    # The current position of fixed
    pos_fixed := function(l) return Position(l, fixed); end;

    # We swap each element of list to the right offset from fixed
    for elem in [1..n] do
        curpos := function(l) return Position(l, elem); end;

        # We should be elem-fixed after fixed
        correct_offset := (elem-fixed) mod Length(l);
        current_offset := function(l) return (curpos(l) - pos_fixed(l)) mod Length(l); end;

        # There are two ways to go, left or right. You can work out
        # which way is faster, but you can also just try both and see
        # which was quicker.

        # First try left.
        left_list := ShallowCopy(list);
        left_moves := 0;
        while current_offset(left_list) <> correct_offset do
            swap_left(left_list, Position(left_list, elem));
            left_moves := left_moves + 1;
        od;

        # Then right.
        right_list := ShallowCopy(list);
        right_moves := 0;
        while current_offset(right_list) <> correct_offset do
            swap_right(right_list, Position(right_list, elem));
            right_moves := right_moves + 1;
        od;

        # Pick whichever was better.
        if left_moves < right_moves then
            list := left_list;
            nswaps := nswaps + left_moves;
        else
            list := right_list;
            nswaps := nswaps + right_moves;
        fi;
    od;

    return nswaps;
end;

SwapsToSort := function(l)
    local n;
    # Just try everything and get the minimum. The point is that one
    # element is already in the right place and we should sort
    # everything around it, but we don't know which one. We have to
    # just try everything. This is a very brute-force stupid way to do
    # this, but it works ok.
    n := Length(l);
    return Minimum(List([1..n], i -> SwapsToSortWithFixed(l, i)));
end;

# Number of adjacent transpositions of elements you need to go from
# perm1 to perm2
AdjacentTranspositionsBetween := function(perm1, perm2)
    local list1, list2, n;

    list1 := MyListPerm(perm1);
    list2 := MyListPerm(perm2);
    n := Length(list1);

    # We rename elements so that list1 "is" [1..n], so we just need to
    # cyclically sort list2.
    list2 := List(list2, elem -> Position(list1, elem));

    # We shift list2 so that 1 is at the start and we don't have to
    # move it. Now we need to count swaps needed to sort list2.
    list2 := ShiftLeft(list2, Position(list2, 1)-1);

    return SwapsToSort(list2);
end;

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

# TODO: make the next bit take irreps from Sage instead of using
# cyclotomics irreps from inside GAP.

# See https://homepages.cwi.nl/~lex/files/symm.pdf for the method we
# now apply to get a smaller semidefinite program.

# First, we block diagonalize action_hom
block_diag_info := BlockDiagonalRepresentationFast(action_hom);
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
nice_basis := Basis(centralizer, norm_cent_basis);

# The multiplication params are defined by B_i B_j = \sum_k \lambda_{i,j}^k B_k
# The convention I use is that mult_param[i][j][k] = \lambda_{i,j}^k
mult_param := NullMat(d, d);
for i in [1..d] do
    for j in [1..d] do
        mult_param[i][j] := Coefficients(nice_basis,
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

# Now use Sage to solve the SDP problem from the paper to get \alpha_m
