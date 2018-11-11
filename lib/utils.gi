# Returns list of first n elems from list
Take@ := function(list, n)
    local result, count;
    result := [];
    count := 0;
    while count < n do
        Add(result, list[count+1]);
        count := count + 1;
    od;
    return result;
end;

# Returns list of all but first n elems from list
Drop@ := function(list, n)
    local result, count, elem;
    result := [];
    count := 0;
    for elem in list do
        if count >= n then
            Add(result, elem);
        fi;
        count := count + 1;
    od;
    return result;
end;

# Returns a list consisting of n copies of elem
Replicate@ := function(elem, n)
    local result, i;
    result := [];
    for i in [1..n] do
        Add(result, elem);
    od;
    return result;
end;

# Takes a list of blocks (possibly different sizes) and constructs a
# block diagonal matrix with those blocks.
InstallGlobalFunction( BlockDiagonalMatrix, function(blocks)
    local combine_blocks, result, block;

    # Combines two blocks into a block diagonal matrix
    combine_blocks := function(b1, b2)
        local len1, len2, new_b1, new_b2;
        len1 := Length(b1);
        len2 := Length(b2);

        # Add len2 zeroes to the end of each row in b1
        new_b1 := List(b1, row -> Concatenation(row, Replicate@(0, len2)));

        # Add len1 zeroes to the start of each row in b2
        new_b2 := List(b2, row -> Concatenation(Replicate@(0, len1), row));

        return Concatenation(new_b1, new_b2);
    end;

    result := [];

    for block in blocks do
        result := combine_blocks(result, block);
    od;

    return result;
end );
