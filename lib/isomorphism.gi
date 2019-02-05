# Various representation isomorphism functions e.g. testing for
# isomorphism and computing an explicit isomorphism.

# Finds the fixed point space of the map A -> tau(gen)*A*rho(gen^-1)
FixedSpace@ := function(rho, tau, gen)

end;

# Picks a random (nonzero) vector in the intersection of some vector
# spaces over C
RandomVectorIntersection@ := function(spaces)

end;

# Wraps an n^2 long list into a n long list of n long lists
WrapMatrix@ := function(vec, n)
    return List([0..n-1], i -> vec{[1+n*i..n*(i+1)]});
end;

InstallGlobalFunction( RepresentationIsomorphism, function(rho, tau)
    local G;

    if not AreRepsIsomorphic(rho, tau) then
        return fail;
    fi;

    G := Source(rho);
    n := DegreeOfRepresentation(rho);

    # We want to find a matrix A s.t. tau(g)*A = A*rho(g) for all g in
    # G. We do this by finding fixed points of the linear maps A ->
    # tau(g)*A*rho(g^-1) (eigenspaces for eigenvalue 1), then we
    # intersect the spaces and pick a random vector (matrix). It will
    # probably be invertible, if not try again.

    # Note: we only need to do this for the generators: homomorphism
    # properties mean this will then work for any g in G
    gens := GeneratorsOfGroup(G);

    fixed_spaces := List(gens, gen -> FixedSpace@(rho, tau, gen));

    repeat
        A_cand := WrapMatrix@(RandomVectorIntersection@(fixed_spaces), n);
    until A_cand^-1 <> fail;

    return A_cand;
end );

# Tells you if two representations of the same group are isomorphic by
# examining characters
InstallGlobalFunction( AreRepsIsomorphic, function(rep1, rep2)
    local G, irr_chars;

    if Source(rep1) <> Source(rep2) then
        return false;
    fi;

    G := Source(rep1);
    irr_chars := IrrWithCorrectOrdering@(G);

    # Writes the characters in the irr_chars basis, they are the same
    # iff they are isomorphic
    return DecomposeCharacter@(rep1, irr_chars) = DecomposeCharacter@(rep2, irr_chars);
end );
