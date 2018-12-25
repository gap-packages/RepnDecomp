Read("utils.g");

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
