InstallGlobalFunction( GroupSumBSGS, function(G, summand)
    local H, n, chain, groups, m, sum, i, cosets, iso, zero, g, S;

    # stab chain computation only works on permutation groups
    iso := IsomorphismPermGroup(G);
    H := Image(iso, G);
    chain := StabChain(H);
    groups := ListStabChain(chain);

    # groups[i] is the stabiliser of the first i-1 points, so
    # groups[1] is H and so on. We start with the last groups[]
    # element, which we have to sum naively, then move to the larger
    # group and so on
    m := Length(groups);

    # summand(g) might be huge to store in memory, recall the whole
    # point of this is to save memory, so only store one image at most
    zero := Zero(summand(One(G)));
    sum := zero;
    for g in GroupStabChain(groups[m-1]) do
        sum := sum + summand(PreImage(iso, g));
    od;

    for i in [m-2,m-3..1] do
        # `sum` is the sum for G_{i+1}

        cosets := RightCosets(GroupStabChain(groups[i]),
                              GroupStabChain(groups[i+1]));

        sum := zero;
        for S in cosets do
            sum := sum + sum * summand(PreImage(iso, Representative(S)));
        od;
    od;

    return sum;
end );
