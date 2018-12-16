# Some algorithms to do with permutations that are needed to compute
# bounds on the crossing number of K_{m, n}

LoadPackage("RepnDecomp");

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

# Converts a permutation to a list. i.e. if you have (a,b,c), the list
# is [a,b,c]. I assume we only get perms of the form (1.....), that is
# cycles with a 1 in them.
MyListPerm := function(perm)
    local result, current;
    result := [];
    current := 1;

    repeat
        Add(result, current);
        current := current^perm;
    until current = 1; # when we reach 1, we have written down the
                       # whole cycle

    return result;
end;

Drop := Drop@RepnDecomp;
Take := Take@RepnDecomp;

# Cyclically moves elements in a list left by n elements
ShiftLeft := function(list, n)
    local shift;
    shift := n mod (Length(list));
    return Concatenation(Drop(list, shift), Take(list, shift));
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

    return Inversions(list2);
end;

# TODO: Something here is incorrect. Make it correct!

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
