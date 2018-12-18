# Some algorithms to do with permutations that are needed to compute
# bounds on the crossing number of K_{m, n}

LoadPackage("RepnDecomp");

Read("compute_q.g");

# Number of inversions in the list. Note that a list is cyclically
# sorted iff it has 0 (sorted) or 1 (cyclically sorted) inversions.
Inversions := function(list)
    local n, image, nswaps, i, j;

    n := Length(list);
    nswaps := 0;

    for i in [1..n] do
        for j in [i+1..n] do
            if list[i] > list[j] then
                nswaps := nswaps + 1;
            fi;
        od;
    od;

    return nswaps;
end;

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
Q := function(p, q) return NumberInterchangesBetween(p, q^(-1)); end;

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
            big[i][j] := NumberInterchangesBetween(ncycles[i], ncycles[j]^-1);
        od;
    od;

    # gives the matrix and the indexing list, so Qmatrix[i][j] =
    # Q(index[i], index[j]).
    return rec(matrix := big,
               index := ncycles);
end;

# Here we fix some number, the end goal is to calculate \alpha_m
m := 5;

# The dimension of the matrices (number of m-cycles)
d := Factorial(m-1);

# This is the group Q is invariant under
G := DirectProduct(SymmetricGroup(m), SymmetricGroup(2));

# a random m cycle
mcycle := MappingPermListList(ShiftLeft([1..m], 1), [1..m]);

# m-cycles index the rows and cols of Q
mcycles := List(mcycle^G);

# This is how G acts on mcycles
action := function(cycle, g)
    local g1, g2, result;

    # Acts via conjugation
    g1 := Image(Projection(G, 1), g);

    # Acts by inverting the cycle
    g2 := Image(Projection(G, 2), g);

    result := cycle;

    # if it's nontrivial, invert
    if g2 <> One(SymmetricGroup(2)) then
        result := result^-1;
    fi;

    result := g1^-1 * result * g1;

    return result;
end;

action_hom := ActionHomomorphism(G, mcycles, action);

# We want the action to be represented as permutation matrices
# Conjugating by any of these matrices fixes Q
action_hom := ConvertRhoIfNeeded@RepnDecomp(action_hom);

# Now we compute the orbits of G on mcycles x mcycles (pi x pi) We do
# this following the paper. Pick some element of pi x pi, keep acting
# on it by the generators g_i until
