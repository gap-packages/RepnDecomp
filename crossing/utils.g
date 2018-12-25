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
