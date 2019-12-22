# Various representation isomorphism functions e.g. testing for
# isomorphism and computing an explicit isomorphism.

# Wraps an n^2 long list into a n long list of n long lists
WrapMatrix@ := function(vec, n)
    local result, current_row, elems_seen;
    result := [];
    current_row := [];
    elems_seen := 0;
    while elems_seen < Length(vec) do
        elems_seen := elems_seen + 1;
        Add(current_row, vec[elems_seen]);
        if Length(current_row) = n then
            Add(result, current_row);
            current_row := [];
        fi;
    od;
    return result;
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
    local G, n, matrix_basis, vector_basis, alpha, triv, fixed_space, A, A_vec, rho_cent_basis, tau_cent_basis, triv_proj, used_tensors, rho_dual, class1, class2, candidate_map, classes, tries, sum, v, v_0, im, orbit, gen, i, rand;

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

    rho_dual := FuncToHom@(G, g -> TransposedMat(Image(rho, g^-1)));

    # The representation alpha : G -> GL(V) (V is the space of
    # matrices).  We only actually need this if we are *not* using
    # Kronecker products.
    alpha := TensorProductOfRepresentations(tau, rho_dual);

    # the projection of V onto V_triv, the trivial canonical summand,
    # is just given by the sum over whole group of alpha(g)

    if ValueOption("use_kronecker") = true then
        # We do this with the BSGS method, this is probably fast. We try
        # to sum in a way that doesn't require us to store any huge
        # Kronecker products.
        triv_proj := GroupSumBSGS(G, g -> KroneckerProduct(Image(tau, g), Image(rho_dual, g)));
    fi;

    classes := ConjugacyClasses(G);

    # The idea here is that we want an element that is fixed by alpha,
    # the natural choice is the sum of an orbit of some random vector
    # v. We know that some choice of v gives an orbit sum that is
    # invertible, so "almost all" choices of v work.

    A := NullMat(n, n);

    tries := 0;

    repeat
        if ValueOption("use_orbit_sum") = true then
            v_0 := RandomInvertibleMat(n);
            sum := v_0;
            orbit := [v_0];

            # This sums the orbit of v_0 under alpha. We know this will
            # terminate since G is finite.
            i := 1;
            while i <= Length(orbit) do
                for gen in GeneratorsOfGroup(G) do
                    im := Image(alpha, gen) * orbit[i];
                    if not (im in orbit) then
                        Add(orbit, im);
                    fi;
                od;
                i := i + 1;
            od;

            A := Sum(orbit);
        elif ValueOption("use_kronecker") = true then
            A := WrapMatrix@(triv_proj * Flat(RandomInvertibleMat(n)), n);
        else
            # Anything involving kronecker products has O(degree^4)
            # space usage, Could be excessive. in some cases we have
            # to resort to summing over the group which usually works
            # but is slow. We can use linearity of alpha and the pair
            # representation of tensor products to avoid excessive
            # memory usage.
            rand := RandomInvertibleMat(n);
            A := Sum(G, g -> Image(alpha, g) * rand);
        fi;
        tries := tries + 1;
    until RankMat(A) = n; # i.e. until A is invertible

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
    until IsLinearRepresentationIsomorphism(candidate, rho, tau);

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
    return IrrVectorOfRepresentation@(rep1, irr_chars) = IrrVectorOfRepresentation@(rep2, irr_chars);
end );

# checks if A rho(g) = tau(g) A
InstallGlobalFunction( IsLinearRepresentationIsomorphism, function(A, rho, tau)
    local gens;
    gens := GeneratorsOfGroup(Source(rho));
    return RankMat(A) = DegreeOfRepresentation(rho) and ForAll(gens, g -> A * Image(rho, g) = Image(tau, g) * A);
end );
