# Various representation isomorphism functions e.g. testing for
# isomorphism and computing an explicit isomorphism.

# Wraps an n^2 long list into a n long list of n long lists
WrapMatrix@ := function(vec, n)
    return List([0..n-1], i -> vec{[1+n*i..n*(i+1)]});
end;

# Gives the full list of the E_ij standard basis matrices for M_n(C)
MatrixBasis@ := function(n)
    local coords, make_mat;
    coords := Cartesian([1..n], [1..n]);

    make_mat := function(coord)
        local ret;
        ret := NullMat(n, n);
        ret[coord[1]][coord[2]] := 1;
        return ret;
    end;

    return List(coords, make_mat);
end;

# Calculates the isomorphism using the cool fact about the product
# (see below)
InstallGlobalFunction( LinearRepresentationIsomorphism, function(rho, tau, args...)
    local G, n, matrix_basis, vector_basis, alpha, triv, fixed_space, A, A_vec, rho_cent_basis, tau_cent_basis, triv_proj, used_tensors, rho_dual;

    # if set, these bases for the centralisers will be used to avoid
    # summing over G
    rho_cent_basis := fail;
    tau_cent_basis := fail;

    if Length(args) >= 2 then
        rho_cent_basis := args[1];
        tau_cent_basis := args[2];
    fi;

    if not AreRepsIsomorphic(rho, tau) then
        return fail;
    fi;

    G := Source(rho);
    n := DegreeOfRepresentation(rho);

    # We want to find a matrix A s.t. tau(g)*A = A*rho(g) for all g in
    # G. We do this by finding fixed points of the linear maps A ->
    # tau(g)*A*rho(g^-1). This is done by considering the
    # representation alpha: g -> (A -> tau(g)*A*rho(g^-1)) and finding a
    # vector in the canonical summand corresponding to the trivial
    # irrep. i.e. a vector which is fixed by all g, which is exactly
    # what we want.

    # The trick we use to avoid summing over G is to notice that alpha
    # is actually \tau \otimes \rho^*, i.e. g -> tau(g) \otimes
    # \rho(g^-1)^T.

    # the representation alpha : G -> GL(V) (V is the space of matrices)
    # TODO: investigate if we can avoid huge matrices in some way

    rho_dual := FuncToHom@(G, g -> TransposedMat(Image(rho, g^-1)));
    alpha := TensorProductOfRepresentations(tau, rho_dual);

    # the projection of V onto V_triv, the trivial canonical summand,
    # is just given by the sum over whole group of alpha(g)
    #
    # we can get a projection into the same space by doing the sum
    # over g in G, h in G of tau(g) \otimes rho^*(h), and we can
    # calculate these group sums using the centraliser bases

    used_tensors := true;

    # if there was no basis, just sum over the whole group
    if rho_cent_basis = fail or tau_cent_basis = fail then
        used_tensors := false;
        alpha := FuncToHom@(G, g -> KroneckerProduct(Image(tau, g), TransposedMat(Image(rho, g^-1))));
        triv_proj := Sum(G, g -> Image(alpha, g));
    else
        # we can just do (sum_{g in G} tau(g)) \otimes (sum_{g in G} rho^*(g))
        # which still gives a projection

        # The group sum for rho^* is the same as for rho, but
        # transposed (the relabelling g -> g^-1 is just a bijection and ^T is linear)
        triv_proj := TensorProductOfMatrices(GroupSumWithCentralizer@(tau, tau_cent_basis),
                                             TransposedMat(GroupSumWithCentralizer@(rho, rho_cent_basis)));

    fi;

    A := NullMat(n, n);

    # Keep picking matrices until we get an invertible one. This would
    # happen with probability 1 if we really picked uniformly random
    # vectors over a ball in C^n^2.
    repeat
        # we pick a "random vector" and project it to get a fixed one
        if used_tensors then
            A := triv_proj * RandomInvertibleMat(n);
        else
            A_vec := triv_proj * Flat(RandomInvertibleMat(n));
            A := WrapMatrix@(A_vec, n);
        fi;
    until RankMat(A) = n;

    return A;
end );

# Calculate the iso without using the cool fact about the conjugation
# action being given by \tau \otimes \rho^*
LinearRepresentationIsomorphismNoProduct@ := function(rho, tau)
    local G, n, matrix_basis, vector_basis, alpha, triv, fixed_space, A, A_vec, alpha_f;

    if not AreRepsIsomorphic(rho, tau) then
        return fail;
    fi;

    G := Source(rho);
    n := DegreeOfRepresentation(rho);

    # We want to find a matrix A s.t. tau(g)*A = A*rho(g) for all g in
    # G. We do this by finding fixed points of the linear maps A ->
    # tau(g)*A*rho(g^-1). This is done by considering the
    # representation alpha: g -> (A -> tau(g)*A*rho(g^-1)) and finding a
    # vector in the canonical summand corresponding to the trivial
    # irrep. i.e. a vector which is fixed by all g, which is exactly
    # what we want.

    # Standard basis for M_n(C)
    matrix_basis := MatrixBasis@(n);

    # The unwrapped vector versions, these are the what the matrices
    # of our representation will act on
    #
    # We construct them in this way to keep track of the
    # correspondence between the matrices and vectors
    vector_basis := List(matrix_basis, Flat);

    # This is the function that gives the representation alpha
    alpha_f := function(g)
        local matrix_imgs, vector_imgs;
        matrix_imgs := List(matrix_basis, A -> Image(tau, g) * A * Image(rho, g^-1));

        # unwrap them and put them in a matrix
        vector_imgs := List(matrix_imgs, Flat);

        # we want the images to be columns not rows (we act from the
        # left), so transpose
        return TransposedMat(vector_imgs);
    end;

    # the representation alpha
    alpha := FuncToHom@(G, alpha_f);

    # trivial representation on G
    triv := FuncToHom@(G, g -> [[1]]);

    # space of vectors fixed under alpha
    fixed_space := IrrepCanonicalSummand@(alpha, triv);

    A := NullMat(n, n);

    # Keep picking matrices until we get an invertible one. This would
    # happen with probability 1 if we really picked uniformly random
    # vectors over C.
    repeat
        # we pick a "random vector"
        A_vec := Sum(Basis(fixed_space), v -> Random(Integers)*v);
        A := WrapMatrix@(A_vec, n);
    until RankMat(A) = n;

    return A;
end;


# checks if it is the case that for all g in G, tau(g)*A = A*rho(g)
IsRepresentationIsomorphism@ := function(rho, tau, A)
    local G, gens, n;

    G := Source(rho);
    n := DegreeOfRepresentation(rho);

    # we only need to check the property for the generators
    gens := GeneratorsOfGroup(G);

    # need bijection and G-action preserving
    return RankMat(A) = n and ForAll(gens, g -> Image(tau, g) * A = A * Image(rho, g));
end;

# calculates an isomorphism between rho and tau by summing over G
# (slow, but works)
# TODO: can I use a trick to sum over the generators instead of G?
InstallGlobalFunction( LinearRepresentationIsomorphismSlow, function(rho, tau, args...)
    local G, n, candidate, tries;
    G := Source(rho);
    n := DegreeOfRepresentation(rho);

    if not AreRepsIsomorphic(rho, tau) then
        return fail;
    fi;

    # we just pick random invertible matrices and sum over the group
    # until we actually get a representation isomorphism. This almost
    # always happens, so we "should" get one almost always on the
    # first time.
    candidate := NullMat(n, n);
    tries := 0;
    repeat
        candidate := RandomInvertibleMat(n);
        candidate := Sum(G, g -> Image(tau, g) * candidate * Image(rho, g^-1));
        tries := tries + 1;
    until IsRepresentationIsomorphism@(rho, tau, candidate);

    if Length(args) > 0 and args[1] = "print tries" then
        Print(tries, " tries\n");
    fi;

    return candidate;
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

# checks if A rho(g) = tau(g) A
InstallGlobalFunction( IsLinearRepresentationIsomorphism, function(A, rho, tau)
    local gens;
    gens := GeneratorsOfGroup(Source(rho));
    return ForAll(gens, g -> A * Image(rho, g) = Image(tau, g) * A);
end );
