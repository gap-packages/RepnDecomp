InstallGlobalFunction( GroupSumBSGS, function(G, summand)
    local H, n, chain, groups, m, sum, i, cosets, iso, zero, g, S, right_reps;

    # stab chain computation only works on permutation groups
    iso := IsomorphismPermGroup(G);
    H := Image(iso, G);
    chain := StabChain(H);
    groups := ListStabChain(chain);

    # groups[i] is the stabiliser of the first i-1 points, so
    # groups[1] is H and so on. We start with the last groups[]
    # element, which we have to sum naively, then move to the larger
    # group and so on
    m := Length(groups)-1;

    groups := List([1..m], i -> GroupStabChain(groups[i]));

    # summand(g) might be huge to store in memory, recall the whole
    # point of this is to save memory, so only store one image at most
    zero := Zero(summand(One(G)));

    # We must sum G_m naively, this should be cheap (it's the
    # stabiliser of a "large" set so should be "small")
    sum := Sum(groups[m], g -> summand(PreImage(iso, g)));

    for i in [m-1,m-2..1] do
        # `sum` is the sum for G_{i+1}

        right_reps := RightTransversal(groups[i], groups[i+1]);
        sum := Sum(right_reps, r -> sum * summand(PreImage(iso, r)));

        # now `sum` is the sum for G_i
    od;

    return sum;
end );
