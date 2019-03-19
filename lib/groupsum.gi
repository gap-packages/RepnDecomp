InstallGlobalFunction( GroupSumBSGS, function(rho)
    local G, H, n, chain, groups, m, sum, i, cosets, iso;

    G := Source(rho);
    n := DegreeOfRepresentation(rho);

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

    sum := Sum(GroupStabChain(groups[m-1]), g -> Image(rho, PreImage(iso, g)));

    for i in [m-2,m-3..1] do
        # `sum` is the sum for G_{i+1}

        cosets := RightCosets(GroupStabChain(groups[i]),
                              GroupStabChain(groups[i+1]));

        sum := Sum(cosets, S -> sum * Image(rho, PreImage(iso, Representative(S))));
    od;

    return sum;
end );
