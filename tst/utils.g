# These are functions for testing, they don't implement anything
# interesting

# Check a subspace is really G-invariant (action given by rho)
IsGInvariant := function(rho, space)
    local G, v, g;

    G := Source(rho);

    if not IsVectorSpace(space) then
        return fail;
    fi;

    for v in Basis(space) do
        for g in GeneratorsOfGroup(G) do
            if not v * Image(rho, g) in space then
                return false;
            fi;
        od;
    od;
    return true;
end;
