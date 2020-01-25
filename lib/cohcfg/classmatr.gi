# BUILD AN IsCohCfg@ OBJECT
CohCfgFromPermGroup@ := function(arg) # computes the number of orbitals, and a transversal.
# arg[1] = G, the group;
# arg[2] (optional) = Omega, the set on which G acts;
# arg[3] (optional) = H, a list of stabilizer subgroups of G, one for each orbit
    # os: (within a loop over I) a list whose x-th sublist is a list of representatives of the H[I]-orbits contained within the x-th G-orbit
    # x: temp int
    # ctr: (within a loop over I) a list showing how many H[I]-orbits each G-orbit contains
    # it: (within a loop) temporary version of i
    # deg: a list of the subdegrees of G on the I-th orbit, that is, the cardinalities of the H[I]-orbits, which correspond to the orbitals of G on the I-th orbit of G acting on Omega (i.e. the I-th transitive component of G acting on Omega)
    # gpairs: a global pairing map, i.e., maps an orbital's global index to that of its transpose
    # ttorbs: keeps track of the total number of orbitals found so far after each row-block of the matrix structure is searched
    # topair:
    # len: (within a loop over I) the number of H[I]-orbits
    # onum: (within a loop over I) a list showing, for each H[I]-orbit, what number it is (in the order it was found) inside its corresponding G-orbit
    # gens: generators of G
    # w: (within a loop) an OrbitNumbers object representing the orbits of the action on the entire set Omega of the stabilizer of the I-th orbit of G
    # onums: a list whose I-th element
    # n: size of Omega, the actioned set
    # G: the underlying group
    # H: list of point stabilizers, one for each orbit
    # orbs: GRAPE-style OrbitNumbers object for G acting on Omega
    # O2
    # d
    # I
    # i
    # J
    # j
    # j0
    # sv: Schreier vector
    local os, x, ctr, it, deg,
        gpairs, ttorbs, topair,
        len, onum, gens, w, onums, G, H, Omega, orbs, O2, d, I, i, J, j, j0, sv;
    G := arg[1]; # set G as the underlying group

    # Sanity check
    if not IsPermGroup(G) then
        Error("CohCfgFromPermGroup@: G must be a permutation group (IsPermGroup object).\n");
    fi;

    if IsBound(arg[2]) then
        if IsHomogeneousList(arg[2]) then
            Omega := arg[2]; # Omega should be a set of points
        elif IsPosInt(arg[2]) then
            Omega := [1..arg[2]]; # but if it is just a number, we assume the range from 1 to the number is intended
        else
            Error("CohCfgFromPermGroup@: Omega must be a set on which G acts (IsHomogeneousList object containing IsPosInt objects) or a natural number (IsPosInt object) n which will be interpreted as [1..n].\n");
        fi;
    else
        Omega := [1..LargestMovedPoint(G)]; # autodetect Omega when it is not given
    fi;

    orbs := OrbitsInfo@(G, Omega); # a GRAPE-style "Numbers"-type record; .orbits is a list of the orbits, .mapping is a list mapping group elements to the number of the orbit they are in, and .representatives is a transversal of the orbits.
    #orbs := List(orbs.representatives, w -> Orbit(G, w));
    d := Length(orbs.representatives); # the size of the transversal set is the number of orbits.

    # check if H[i] (the stabilizers of representatives of orbits) were given, and if so, make sure they are appropriate
    # XXX: is the orbit ordering really going to be the same every time the orbits are generated?
    if IsBound(arg[3]) then
        if IsList(arg[3]) then
            H := arg[3]; # H should be a list, but if it's not, assume it's a singleton list (transitive group?)
        else
            H := [arg[3]];
        fi;

        if Length(H) <> d then # check if the number of point stabilizers matches the number of orbits calculated
            Error("CohCfgFromPermGroup@: The wrong number of point stabilisers was supplied.\n");
        fi;

        for i in [1..d] do # for each orbit...
            if ForAny(GeneratorsOfGroup(H[i]), x -> orbs.representatives[i]^x <> orbs.representatives[i]) then # check that the point stabilizer given is truly a stabilizer for the points in that orbit
                Error("CohCfgFromPermGroup@: H[", i, "] does not fix the point ", orbs.representatives[i], " .\n");
            fi;
        od;
    else
        H := List(orbs.representatives, i -> Stabilizer(G, i)); # otherwise, just grow the stabilizers by finding the stabilizers for each representative in the transversal set of the orbits
    fi;

    onums := [];
    ttorbs := 0;

    for I in [1..d] do # for each I-th orbit...
        #Add(onums, List([1..d], rec(orbnums := OrbitNumbers(H[i], n), paired := [])));
        w := OrbitsInfo@(H[I], Omega); # w is an OrbitNumbers object representing the orbits of the action on the entire set Omega of the stabilizer of the I-th orbit of G

        # need to localize orbit representatives per G-orbit.
        ctr := List([1..d], i -> 0); # will be a list showing how many H[I]-orbits each G-orbit contains
        len := Length(w.representatives); # len is the number of H[I]-orbits
        onum := EmptyPlist(len); # will be a list showing, for each H[I]-orbit, what number it is (in the order it was found) inside its corresponding G-orbit
        deg := List([1..len], i -> 0); # deg will become a list of the subdegrees of G on the I-th orbit, that is, the cardinalities of the H[I]-orbits, which correspond to the orbitals of G on the I-th orbit of G acting on Omega (i.e. the I-th transitive component of G acting on Omega)

        # computing subdegrees by counting the elements in each H[I]-orbit
        for i in w.mapping do
            deg[i] := deg[i] + 1;
        od;

        os := List([1..d], i -> []); # will be a list whose x-th sublist is a list of representatives of the H[I]-orbits contained within the x-th G-orbit
        for i in [1..len] do # for each i-th H[I]-orbit...
            x := orbs.mapping[w.representatives[i]]; # find which x-th G-orbit of Omega the i-th H[I]-orbit is a subset of
            Add(os[x], w.representatives[i]);          # add the representative of the i-th H[I]-orbit to the x-th sublist of os
            ctr[x] := ctr[x] + 1;
            onum[i] := ctr[x];
        od;
        # now, os is a list whose x-th sublist is a list of representatives of the H[I]-orbits contained within the x-th G-orbit
        # now, ctr is a list showing how many H[I]-orbits each G-orbit contains
        # now, onum is a list showing, for each H[I]-orbit, what number it is (in the order it was found) inside its corresponding G-orbit

        ### internally, orbitals are indexed by pairs (G-orbit, H[I]-orbit)
        # internally, orbitals are indexed by pairs (I, J) where I is a G-orbit and J is an H[I]-orbit of Omega
        Add(onums, rec(
            globalTwoOrbitNumbers := List(w.mapping, i -> i + ttorbs), # a list mapping elements alpha in Omega to the global orbital number of the orbital containing (rep, alpha), where rep is the representative of the I-th G-orbit
            stabiliser := w, # an OrbitNumbers object representing the orbits of the action on the entire set Omega of the stabilizer of the I-th orbit of G
            representatives := os, # a list whose x-th sublist is a list of representatives of the H[I]-orbits contained within the x-th G-orbit
            perBlockOrbitNumber := onum,
            numberTwoOrbits := Length(Flat(os)), # number of orbitals found rooted in the I-th G-orbit (i.e., H[I]-orbits of Omega)
            paired := List([1..d], x -> []), # pairs the orbital matrices with their transposes
            subdegree := deg
        )); # nonzero row sums of orbital matrices
        ttorbs := ttorbs + len; # increment the uid offset for counting H[?]-orbits
    od;
    # now, onums is a list of records, per G-orbit on Omega, with the I-th list containing:
    #   .globalTwoOrbitNumbers - gives an OrbitsInfo@.mapping for H[I] acting on Omega, but with the indices of the (H[I]-)orbits shifted to be unique over all I
    #   .stabiliser - an OrbitsInfo@ object representing the H[I]-orbits on Omega
    #   .representatives - a list whose x-th sublist is a list of representatives of the H[I]-orbits contained within the x-th G-orbit
    #   .perBlockOrbitNumber - a list showing, for each H[I]-orbit, what number it is (in the order it was found) inside its corresponding G-orbit
    #   .numberTwoOrbits - the total number of H[I]-orbits (which is equal to the total number of orbitals of G whose first coordinate is restricted to the I-th G-orbit)
    #   .paired -
    #   .subdegree - a list of the subdegrees of G acting on the I-th G-orbit on Omega
    # now, ttorbs is the total number of H[I] orbits over all I. (which is equal to the total number of orbitals of G)
        if Omega<>[1..Length(Omega)] then
            Error("CohCfgFromPermGroup@: Schreier Vector of a non-solid range is not implemented.\n");
        fi;
    sv := NullGraph(G,Length(Omega)).schreierVector; # the Schreier vector for G, computed by GRAPE's NullGraph function
    gens := GeneratorsOfGroup(G);

    # establishing pairings
    gpairs := List([1..ttorbs], i -> 0); # will be a global pairing map, i.e., maps an orbital's global index to that of its transpose

    for I in [1..d] do # for each I-th G-orbit...
        i := orbs.representatives[I]; # the representative of the orbit
        for J in [1..d] do # for each J-th G-orbit...
            for j0 in [1..Length(onums[I].representatives[J])] do
                j := onums[I].representatives[J][j0]; ### global column number
                topair := onums[I].globalTwoOrbitNumbers[j];
                ### need to canonise (j, i), i.e. find g in G that maps j to orbs.representatives[J].
                ### (using the Schreier vector sv in the G's nullgraph for this)
                ### then g(i) should be located in the orbits of the stabilizer of
                ### orbs.representatives[J], i.e. in onums[J][I].stabiliser.mapping .
                ### So onums[I].paired[J][j0] must get the corresponding entry assigned.

                w := sv[j];
                it := i;
                while w > 0 do
                    j := j/gens[w];
                    it := it/gens[w];
                    w := sv[j];
                od;

                gpairs[topair] := onums[J].globalTwoOrbitNumbers[it];

                onums[I].paired[J][j0] := onums[J].perBlockOrbitNumber[onums[J].stabiliser.mapping[it]];
                ###Error("2");

                if J > I then ### ??? or onums[J].stabiliser.mapping[i] <> j0 then
                    onums[J].paired[I][
                        onums[J].perBlockOrbitNumber[onums[J].stabiliser.mapping[it]]
                    ] := onums[J].stabiliser.mapping[onums[J].representatives[I][j0]];
                fi;
            od;
        od;
    od;

    return rec(
        isTwoOrbitNumbers := true, # object type indicator
        group := G, # underlying group of the coherent configuration (can this be generalized away?)
        Omega := Omega, # the set Omega, which is acted upon by G
        degree := Length(Set(Omega)), # the degree of the coherent configuration (is this needed?)
        orbitMapping := orbs.mapping, # map list showing which orbit each point is in
            orbitNumbers := ~.orbitMapping, # for backwards-compatibility
        P := [], # to be populated later, will be a 3D array of p_{i,j}^k, the intersection numbers
        orbitNumberRepresentatives := orbs.representatives, # list of representatives of the orbits
        orbitLengths := List(orbs.orbits, Length), # list of lengths of the G-orbits
        twoOrbitNumbers := onums, # a collection of stuff... see above
        schreierVector := sv, # a Schreier vector of the graph
        dimension := ttorbs, # the number of orbitals, which is the same as cardinality of the coherent configuration and the dimension of the algebra it generates
        pairing := PermList(gpairs) # a permutation which, when applied to the ordered list of orbitals, will transpose each orbital
    );
end;

###########################################################################
# standard function to test whether this is the consistent record
###########################################################################
IsCohCfg@ := function(T)
    return IsRecord(T) and IsBound(T.isTwoOrbitNumbers) and T.isTwoOrbitNumbers;
end;

###########################################################################
# underlying group
###########################################################################
CCGroup@ := function(T)
    if not IsCohCfg@(T) then
        Error("CCGroup@: this is not a coherent configuration.\n");
    fi;

    return T.group;
end;

###########################################################################
# size of the domain operated on by the underlying group
###########################################################################
CCDegree@ := function(T)
    if not IsCohCfg@(T) then
        Error("CCDegree@: this is not a coherent configuration.\n");
    fi;

    return T.degree;
end;

###########################################################################
# orbit lengths of the underlying group
###########################################################################
CCFiberSizes@ := function(T)
    if not IsCohCfg@(T) then
        Error("CCFiberSizes@: this is not a coherent configuration.\n");
    fi;

    return T.orbitLengths;
end;

CCCountByLeftFiber@ := function(TN, I) #orbitals on an orbit
    if not IsCohCfg@(TN) then
        Error("CCCountByLeftFiber@: this is not a coherent configuration.\n");
    fi;

    return TN.twoOrbitNumbers[I].numberTwoOrbits;
end;

indexInLeftFiber@ := function(TN, onum, blk)
    local b;

    # TODO: fail if given CC element is in the wrong block

    for b in [1..blk-1] do
        onum := onum - Length(TN.twoOrbitNumbers[b].stabiliser.representatives);
    od;

    return onum;
end;

# convert the rec(orbitNumber, twoOrbitNumber) into the global number
CCGlobalIndex@ := function(TN, lnum)
    local j, s, I;

    if not IsCohCfg@(TN) then
        Error("CCGlobalIndex@: this is not a coherent configuration.\n");
    fi;

    I := lnum.orbitNumber;
    s := lnum.twoOrbitNumber;
    for j in [1..I-1] do
        s := s + CCCountByLeftFiber@(TN, j);
    od;

    return s;
end;

# and the reverse conversion
CCLocalIndex@ := function(TN, i)
    local I, lnum;

    if not IsCohCfg@(TN) then
        Error("CCLocalIndex@: this is not a coherent configuration.\n");
    fi;

    I := 1;
    while i > CCCountByLeftFiber@(TN, I) do
        i := i - CCCountByLeftFiber@(TN, I);
        I := I + 1;
    od;

    lnum := rec(orbitNumber := I, twoOrbitNumber := i);
    return lnum;
end;

########################################################################
# for a 2-orbit t, compute Trace(tt^*)
########################################################################
CCElementSize@ := function(T, t)
    local lon, orblen, deg;

    if not IsCohCfg@(T) then
        Error("CCElementSize@: this is not a coherent configuration.\n");
    fi;

    lon := CCLocalIndex@(T, t);
    orblen := CCFiberSizes@(T)[lon.orbitNumber];
    deg := T.twoOrbitNumbers[lon.orbitNumber].subdegree[lon.twoOrbitNumber];

    return deg*orblen;
end;

########################################################################
# for an arc (a, b), find the (global) number of the two-orbit it lies in
########################################################################
CCElementContainingPair@ := function(TN, pair)
    local sv, w, gens, I, lnum, globalNumber,
        a, b;
    a := pair[1]; b := pair[2];
    if (not IsCohCfg@(TN)
        or a>CCDegree@(TN)
        or b>CCDegree@(TN)
        or not IsPosInt(a)
        or not IsPosInt(b)
    ) then
        Error("CCElementContainingPair@: Usage is CCElementContainingPair@(IsCohCfg@ object, pair := [int, int]).\n");
    fi;

    sv := TN.schreierVector;
    gens := GeneratorsOfGroup(CCGroup@(TN));
    w := sv[a];
    while w > 0 do
        a := a/gens[w];
        b := b/gens[w];
        w := sv[a];
    od;

    I := TN.orbitMapping[a];
    return TN.twoOrbitNumbers[I].globalTwoOrbitNumbers[b];
end;

CCNumFibers@ := function(TN)
    if not IsCohCfg@(TN) then
        Error("CCNumFibers@: Usage is CCNumFibers@(IsCohCfg@ object).\n");
    fi;

    return Length(TN.orbitNumberRepresentatives);
end;

########################################################################
# dimension of the algebra
########################################################################
CCDimension@ := function(TN)
    local x;
    if not IsCohCfg@(TN) then
        Error("CCDimension@: Usage is CCDimension@(IsCohCfg@ object).\n");
    fi;

    return TN.dimension;
end;

#########################################################################
# list a pair for each 2-orbit
#########################################################################
CCTransversal@ := function(T)
    local i, x;

    if not IsCohCfg@(T) then
        Error("CCTransversal@: this is not a coherent configuration.\n");
    fi;

    return Concatenation(List(
        [1..Length(T.orbitNumberRepresentatives)],
        i -> List(
            T.twoOrbitNumbers[i].stabiliser.representatives,
            x -> [T.orbitNumberRepresentatives[i], x]
        )
    ));
end;

#########################################################################
# a naive implementation ignoring block structures
# (initialization (?) )
# 6.12.08 - added sparsity handling
#########################################################################
doCoeffsComputation@ := function(T, ci, cj, I, J, k)
    local t, p, i, M;

    M := List([1..T.dimension], i -> []);
    for i in [1..CCDegree@(T)] do
        p := PositionSet(List(M[ci[i]], t -> t[1]), cj[i]);
        if p = fail then
            AddSet(M[ci[i]], [cj[i], 1]);
        else
            M[ci[i]][p][2] := M[ci[i]][p][2] + 1;
        fi;
    od;

    #M := NullMat(T.dimension, T.dimension);
    #for i in [1..CCDegree@(T)] do
    #   M[ci[i]][cj[i]] := M[ci[i]][cj[i]] + 1;
    #od;

    T.P[k] := M;

    return 1;
end;

#########################################################################
# for a given k, compute p^k_{ij} for all i, j
#########################################################################
doCoeffsPerCCElement@ := function(TN, k)
    local word, w, sv, gens, x, iii, ci, I, twoonum, Ir,
        gamma, cj, J, M, i, j;

    Ir := CCLocalIndex@(TN, k);
    I := Ir.orbitNumber;
    twoonum := Ir.twoOrbitNumber;

    # locate the arc (i, j) from the (I, J)-block
    j := TN.twoOrbitNumbers[I].stabiliser.representatives[twoonum];

    # the k-th 2-orbit is located in (I, J)-block
    J := TN.orbitMapping[j];
    i := TN.orbitNumberRepresentatives[I];

    # now we have to compute the "scalar product" of i-th row and j-th column
    # of the "big matrix"
    # and store it in M. More precisely, M[a][b] will get the number of points
    # p such that the arc (i, p) has colour a and the arc (p, j) - colour b.
    #
    # the i-th row is already there, it is
    # ci=TN.twoOrbitNumbers[I].globalTwoOrbitNumbers;
    #
    # the j-th column need to be computed from
    # TN.twoOrbitNumbers[J].globalTwoOrbitNumbers by applying the suitable
    # permutation and then the pairing.
    ci := TN.twoOrbitNumbers[I].globalTwoOrbitNumbers;

    sv := TN.schreierVector;
    gens := GeneratorsOfGroup(CCGroup@(TN));
    w := sv[j];

    # need to compute the word in generators,
    word := [];
    while w>0 do
        j := j/gens[w];
        Add(word, w);
        w := sv[j];
    od;

    # and apply it to cj
    cj := List(TN.twoOrbitNumbers[J].globalTwoOrbitNumbers, w -> w^TN.pairing);
    for w in Reversed(word) do
        cj := Permuted(cj, gens[w]);
    od;

    doCoeffsComputation@(TN, ci, cj, I, J, k);
    return 1;
end;

CCPopulateCoeffs@ := function(TN)
    local k;

    if IsEmpty(TN.P) then # compute them
        List([1..CCDimension@(TN)], k ->
            doCoeffsPerCCElement@(TN, k)
        );

        return 1;
    fi;

    return 2; # they were already computed
end;

##########################################################################
# retrieving p^a_{bc}
##########################################################################
CCCoeff@ := function(TN, a, b, c) # XXX Make this actually check whether the coefficients have been generated
    local p, x;
    #return TN.P[a][b][c]; # naive, non-sparse implementation
    p := PositionSet(List(TN.P[a][b], x -> x[1]), c);

    if p = fail then
        return 0;
    else
        return TN.P[a][b][p][2];
    fi;
end;

###########################################################################
# computing the isomorphism to the regular representation, i.e.
# A_b -> (L_b)_{ac} := p^a_{bc}, for a, b, c=1, ..., algebra dimension
###########################################################################
CCIntersectionMat@ := function(TN, b)
    local L, d, a, c;

    d := CCDimension@(TN);
    L := NullMat(d, d);
    for a in [1..d] do
        for c in [1..d] do
            L[a][c] := CCCoeff@(TN, a, b, c);
        od;
    od;

    return L;
end;

# retrieve the regular representation matrices
CCIntersectionMats@ := function(TN)
    local b;
    if IsEmpty(TN.P) then
        CCPopulateCoeffs@(TN);
    fi;

    return List([1..CCDimension@(TN)],
        b -> CCIntersectionMat@(TN, b)
    );
end;

# sparse version
CCIntersectionMatLIL@ := function(TN, b)
    local v, a, c, d, L;

    d := CCDimension@(TN);
    L := List([1..d], i -> []);
    for a in [1..d] do
        for c in [1..d] do
            v := CCCoeff@(TN, a, b, c);
            if v <> 0 then
                AddSet(L[a], [c, v]);
            fi;
        od;
    od;

    return L;
end;

# packed version
# for each non-zero entry M_{ij}, list the tuple [i,j,M_{ij}]
# and return the resulting list (ordered 1st by i, and then by j)
CCIntersectionMatCOO@ := function(TN, b)
    local v, a, c, d, L;

    d := CCDimension@(TN);
    L := [];
    for a in [1..d] do
        for c in [1..d] do
            v := CCCoeff@(TN, a, b, c);
            if v <> 0 then
                Add(L,[a,c,v]);
            fi;
        od;
    od;

    return L;
end;

# packed version
CCIntersectionMatsCOO@ := function(TN)
    if IsEmpty(TN.P) then
        CCPopulateCoeffs@(TN);
    fi;

    return List([1..CCDimension@(TN)], x-> CCIntersectionMatCOO@(TN,x));
end;

# transposes a sparse matrix
TransposedMatLILMutable@ := function(A)
    local p, d, i, j, L;

    d := Length(A);
    L := List([1..d], i -> []);
    for i in [1..d] do
        for p in A[i] do
            AddSet(L[p[1]], [i, p[2]]);
        od;
    od;

    return L;
end;

###############################################################
# this is mainly for debugging purposes
###############################################################
CCBasisMats@ := function(TN)
    local M, i, j, n, d;

    if not IsCohCfg@(TN) then
        Error("CCBasisMats@: Usage is CCBasisMats@(IsCohCfg@ object).\n");
    fi;

    if IsBound(TN.basisMats) then
        return TN.basisMats;
    fi;

    d := CCDimension@(TN);
    n := CCDegree@(TN);
    M := List([1..d], x -> NullMat(n, n));
    for i in [1..n] do
        for j in [1..n] do
            M[CCElementContainingPair@(TN, [i, j])][i][j] := 1; # TODO: can this be done more efficiently?
        od;
    od;

    for i in [1..d] do
        if M[i]<>TransposedMat(M[i^TN.pairing]) then
            Error("CCBasisMats@: Stored pairing is faulty.\n");
        fi;
    od;

    MakeImmutable(M);
    TN.basisMats := M;

    return TN.basisMats;
end;

# produce the aggregate basis matrix whose (i,j)-th entry gives the number of the CC element which (i,j) lies in (0-based)
CCAggregateBasisMat@ := function(TN)
    local As, A;

    if not IsCohCfg@(TN) then
        Error("CCAggregateBasisMat@: Usage is CCAggregateBasisMat@(IsCohCfg@ object).\n");
    fi;

    if IsBound(TN.aggregateBasisMat) then
        return TN.aggregateBasisMat;
    fi;

    # TODO: do this more efficiently
    As := CCBasisMats@(TN);
    A := Sum([2 .. CCDimension@(TN)], x -> (x - 1) * As[x]);

    MakeImmutable(A);
    TN.aggregateBasisMat := A;

    return TN.aggregateBasisMat;
end;
