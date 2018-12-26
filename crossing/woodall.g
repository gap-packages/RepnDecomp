# This is an implementation of an algorithm due to D.R. Woodall to
# compute the number of adjacent interchanges needed to change a cycle
# a into another cycle b.

LoadPackage("RepnDecomp");
Read("utils.g");

# TODO: This algorithm is incorrect for some reason, not sure why
NumberInterchangesBetween := function(a, b)
    local al, bl, n, f, dist, alpha, beta, gamma;

    # Convert them to lists
    al := MyListPerm(a);
    bl := MyListPerm(b);

    n := Length(al);

    # We rename elements so that al "is" [0..n-1]. This is so we match
    # up with Woodall's description and avoid off-by-one errors.
    bl := List(bl, elem -> Position(al, elem)-1);
    al := [0..(n-1)];


    # This is a useful function defined by Woodall
    f := function(r)
        if r = 0 then
            return 0;
        elif 0 < r and r < n/2 then
            return 2*r - 1;
        elif r = n/2 then
            return 2*r - 2;
        fi;
        Error("Bad value for r: ", r);
    end;

    # We imagine the entries of a spaced evenly around a circle. There
    # are n different ways to superimpose b on top of a, dist_i
    # calculates the distance assuming one of the n
    # superimpositions. We consider the ith choice to be b, but
    # shifted left i times.

    # alpha(i, j) denotes the number in a coinciding with j in b
    # when b is in position i.
    alpha := function(i, j)
        local posj, pos_shifted;

        # The 0-indexed position of j in b
        posj := Position(bl, j)-1;

        # Now shift left i
        pos_shifted := (posj - i) mod n;

        # The value of a there is known, it's just the position.
        return pos_shifted;
    end;

    # Now j - alpha(i, j) is equal (mod n) to the distance that j
    # has to move in the positive direction from its position in b
    # to its final position

    # beta(i, j) denotes the absolute value of the integer of
    # minimum absolute value that is congruent (mod n) to j -
    # alpha(i, j), so beta(i ,j) is the shortest distance that j
    # must move.
    beta := function(i, j)
        local x;
        x := (j - alpha(i, j)) mod n;
        return Minimum([x, n - x]);
    end;

    gamma := function(i, j) return (1/2)*f(beta(i, j)); end;

    dist := i -> Sum([0..(n-1)], j -> gamma(i, j));

    # The final result is the smallest distance between a and b over
    # all choices of rotation of b.
    return Minimum(List([0..n-1], i -> dist(i)));
end;
