####################################################################
# let g in G, and A_i's the orbitals of G. Computing
# g^G as a linear combination of A_i's.
####################################################################
# ClassSum@(group, conjugacyclass)
# ClassSum@(group, element, sizeofelementsconjugacyclass)
ClassSum@ := function(arg)
    local arcs, oinf, i, g, G, n, len, v, x, y, fps, onums, sonums, m, orbs, T;
    
    T := arg[1];
    G := CCGroup@(T);
    if not IsPerm(arg[2]) then 
        if G <> ActingDomain(arg[2]) then 
            Error("ClassSum@: The second argument must be a conjugacy class of the underlying group of the first argument.\n");
        fi;
        
        g := Representative(arg[2]);
        len := Size(arg[2]);
    else # no checks are made; g is a representative, len the size of g^G
        g := arg[2];
        len := arg[3];
    fi;
    
    n := CCDegree@(T);
    v := 0*[1..T.dimension]; # coefficients of the decomposition
    
#   if 1=0 then # skip this!!!
#       fps := Filtered([1..n], x -> x^g = x); # fixed points 
#       # (essentially, character computation for each G-orbit)
#       onums := List(fps, x -> T.orbitMapping[x]);
#       sonums := Set(onums);
#       m := List(sonums, x -> Size(Filtered(onums, y -> y = x)));
#       orbs := List(sonums, x -> CCElementContainingPair@(T, [T.orbitNumberRepresentatives[x],
#           T.orbitNumberRepresentatives[x]
#       ]));
#
#       # g hits i-th diagonal 2-orbit orbs[i] m[i] times;
#       # there are len elements in g^G, so we get len*m[i] hits;
#       # by symmetry, the contribution of i-th diagonal element is len*m[i]/(length of i-th orbit)
#       
#       for i in [1..Length(m)] do 
#           v[orbs[i]] := len*m[i]/CCFiberSizes@(T)[i];
#       od;
#   fi; # end of "skip this!!!"
    
    # cycles of the permutation
    arcs := List(Cycles(g, [1..n]), x -> rec(
        #o := CCElementContainingPair@(T, [x[1], x[2]]),
        o := CCElementContainingPair@(T, [x[Length(x)], x[1]]),
        l := Length(x))
    );
    
    # g hits i-th 2-orbit arcs[i].o by a cycle of length arcs[i].l
    # there are len elements in g^G, so we get len hits by a arcs[i].l-cycle;
    # by symmetry, the contribution of this 2-orbit is len*arcs[i].l/(# of 1s in its nxn-matrix)
    
    for i in [1..Length(arcs)] do 
        oinf := CCLocalIndex@(T, arcs[i].o);
        x := oinf.orbitNumber;
        v[arcs[i].o] := v[arcs[i].o] + len*arcs[i].l / (
            CCFiberSizes@(T)[x] * T.twoOrbitNumbers[x].subdegree[oinf.twoOrbitNumber]
        );
    od;
    
    return v;
end;

####################################################################
# computing all the classsums and storing them in the corr. record
####################################################################
ClassSums@ := function(arg)
    local T, x, sums, g_i, l_i, cc;
    
    T:=arg[1];
    if not IsBound(T.classSums) then 
        if not IsBound(arg[2]) then 
            cc := ConjugacyClasses(CCGroup@(T));
            g_i := List(cc, Representative);
            l_i := List(cc, Size);
        else 
            g_i := arg[2];
            l_i := arg[3];
        fi;
        
        sums := List([1 .. Length(g_i)], x -> ClassSum@(T, g_i[x], l_i[x]));
        T.classSums := rec(sums := sums, reps := g_i, lengths := l_i);
    fi;
    
    return T.classSums;
end;

####################################################################
# Computing the projection of the regular representation
# onto a representation with the character \chi.
# 
# We want to find
#     (dim(\chi) / |G|) * \sum_C \chi(C)^* \bar C
# Here the index C in "\sum_C" ranges over conjugacy classes of G,
# and \bar C is the sum of the elements in C, after they have been
# taken through the representation from G to the basis algebra of
# the coherent configuration and from there through the isomorphism
# from the basis algebra to the intersection algebra of the CC.
# 
# Note that
#     \sum_C \chi(C)^* \bar C
#   = \sum_C (\chi(C)^* \sum_{i=1}^d ClassSum@(C)_i P_i)
#   = \sum_{i=1}^d (\sum_C ClassSum@(C)_i \chi(C)^*) P_i)
####################################################################
ProjComp@ := function(T, chi)
    local i, v, M, d, s_r, sums, s, c;
    
    s_r := ClassSums@(T);
    
    return (chi[1] / Size(CCGroup@(T))) * Sum([1 .. CCDimension@(T)], i -> 
        Sum([1 .. Length(s_r.reps)], C -> 
            s_r.sums[C][i] * ComplexConjugate(chi[C])
        ) * CCIntersectionMat@(T, i)
    );
end;

####################################################################
# naive diagonalisation (block diagonalization of the regular representation)
####################################################################
ProjIrr@ := function(arg)
    local ct, T, x, p, G, irrs, projs, linindeprows;
    
    linindeprows := function(V)
        local M, r, len, t;
        M := [];
        r := RankMat(V);
        len := 0;
        while len < r do
            t := PositionProperty(V, y -> RankMat(Concatenation(M, [y])) > len);
            Append(M, [V[t]]);
            V := V{[t+1..Length(V)]};
            len := len+1;
        od;
        
        return M;
    end;
    
    T := arg[1];
    CCPopulateCoeffs@(T);
    ClassSums@(T);
    ct := [];
    
    if not IsBound(arg[2]) then 
        G := CCGroup@(T);
        ct := CharacterTable(G);
        irrs := Irr(ct);
    else 
        irrs := arg[2]; # TODO instead of accepting irreducibles as a list, accept a character table and connect it to the group
    fi;
    
#   G := CCGroup@(T);
#   if not IsBound(arg[2]) then 
#       ct := CharacterTable(G);
#   else 
#       ct := arg[2];
#       ConnectGroupAndCharacterTable(G, ct);
#   fi;
#   irrs := Irr(ct);
    
    projs := List(irrs, x -> rec(p := ProjComp@(T, x), dim := x[1]));
    projs := Filtered(projs, x -> x.p <> 0*projs[1].p);
    projs := List(projs, x -> rec(p := x.p, rank := RankMat(x.p), dim := x.dim));
    if Sum(List(projs, G -> G.rank)) <> T.dimension then 
        Error("ProjIrr@: wrong projectors? ");
    fi;
    
    # "stack" the individual X (which trim matrix to its block component)
    return rec( # TODO this record is ugly, make it sensible (xref BlocksOfMat@)
        projmat := Concatenation(List(projs, G -> linindeprows(G.p))),
        blocksizes := List(projs, G -> G.rank),
        dims := List(projs, G -> G.dim),
        ct := ct
    );
end;

####################################################################
# obtaining blocks of a matrix M in the regular representation
####################################################################
BlocksOfMat@ := function(P, M) # TODO coordinate with format returned by ProjIrr@
    local u, v, b, X, i, bsizes, p, V, k;
    
    p := P.projmat*M*P.projmat^-1;
    bsizes := P.blocksizes;
    V := [];
    i := 0;
    k := 0;
    
    for b in bsizes do 
        k := k + 1;
        if b > 0 then 
            X := NullMat(b,b);
            for u in [1..b] do 
                for v in [1..b] do 
                    X[u][v] := p[i + u][i + v];
                od;
            od;
            
            Add(V, rec(blksz := X, char := P.dims[k]));
            i:=i + b;
        fi;
    od;
    
    return V;
end;
