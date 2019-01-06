# These are functions for testing, they don't implement anything
# interesting

LoadPackage("RepnDecomp");

# Check a subspace is really G-invariant (action given by rho)
IsGInvariant := function(rho, space)
    local G, v, g;

    G := Source(rho);

    if not IsVectorSpace(space) then
        return fail;
    fi;

    for v in Basis(space) do
        for g in GeneratorsOfGroup(G) do
            if not Image(rho, g) * v in space then
                return false;
            fi;
        od;
    od;
    return true;
end;

DirectSumRepList := function(reps)
    local G, gens, imgs, H;
    G := Source(reps[1]);
    gens := GeneratorsOfGroup(G);
    imgs := List(gens, g -> BlockDiagonalMatrix(List(reps, rep -> Image(rep, g))));
    H := Group(imgs);
    return GroupHomomorphismByImages(G, H, gens, imgs);
end;
