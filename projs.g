projs := function(irr,rep,a,b) # p_{ab}, after Serre
    local p, g, G, n, W;

    W := Range(irr);  # the irreducible
    n := Length(W.1); # dimension of W
    G := Source(irr);  # the group
    if G <> Source(rep) then 
       Error("projs: groups acting are not the same!");
    fi;
    p := 0*W.1;
    for g in Elements(G) do
       p := p + (n/Order(G))*Image(irr,g^-1)[b][a]*Image(rep,g);
    od;
    return p;
end;


regrep :=function(G, nc) # for testing purposes: nc copies of the regular rep of G
    local r, a, m, x; 
    r := Range(RegularActionHomomorphism(G));
    m := List(GeneratorsOfGroup(r), a->PermutationMat(a, Order(G), Rationals));
    if nc>1 then 
        m := List(m, x-> BlockMatrix(List([1..nc], a-> [a,a,x]), nc, nc));
    fi;
    return GroupHomomorphismByImages(G, Group(m), GeneratorsOfGroup(G), m);
end;


Virr := function(irr, rep) # compute the list of bases of W(x_1^{k})
    local F, p, G, n, W, p_a1, r, B, V1, v;

    W := Range(irr);
    n := Length(W.1);
    F := FieldOfMatrixGroup(Range(rep));
    p_a1 := List([1..n], p->projs(irr, rep, p, 1));
    V1 := VectorSpace(F, p_a1[1]);
    B := [];
    for v in Basis(V1) do
       r := []; 
       for p in p_a1 do
          Add(r, p*v);
       od;
       Append(B, r);
    od;
    return B;
end;


conjmat := function(irrs, rep) # do B*gen*B^-1 to get block-diagonalisation
    local pat, p, G, n, B, irr, d, x, M;

    G := Source(rep);  # the group
    B := [];
    pat :=[]; 
    for irr in irrs do
      M := Virr(irr, rep);
      d := Length(Image(irr).1);
      Append(B, M); 
      Add(pat, rec(mul:=Length(M)/d, dim:=d));
    od;
    return [B,pat];
end; 

sparsify := function(p, M)
    local n, k, B;
    if p=[] then
       return [];
    fi;
    k:=1;
    B:=[];
    for n in p do
        Add(B,M{[k..k+n-1]}{[k..k+n-1]});
        k := k+n;
    od;
    return B;
end;

exppat := function(recs)
    local x, l, r;
    l := [];
    for r in recs do
       Append(l,List([1..r.mul], x->r.dim));
    od;
    return l;
end; 

tomat := function(r, u) # convert a vector u of length r.dim*r.mul to (r.mul x r.dim)-matrix 
    local k, i, M;
    k := 1;
    M:=[];
    for i in [1..r.mul] do
      Add(M,u{[k..k+r.dim-1]});
      k := k+r.dim;
    od;
    return M;
end;

conjrep := function(irrs, rep) # block-diagonal form
    local B, x, G;
    B := conjmat(irrs, rep);
    G := Image(rep);
    return List(GeneratorsOfGroup(G), x-> sparsify(exppat(B[2]), B[1]*x*B[1]^-1));
end;

symmsquare := function(rep)
    local mats, mat, newmat, row, G, gens, dim, nmats, i, j, k, m, n, x;
    mats := GeneratorsOfGroup(Range(rep));
    nmats := Length(mats); 
    dim := Length(mats[1]);
    gens:=[];
    for mat in mats do
      newmat:=[];
      for i in [1..dim] do
         for j in [i..dim] do
            row:=[];
            for k in [1..dim] do
               for m in [k..dim] do
                  x:=mat[m][j] * mat[k][i]+mat[m][i]*mat[k][j];
                  if k<m then 
                      Add (row, x);
                  else
                      Add (row, (1/2)*x);
                  fi;
               od;
            od;
            Add (newmat, row);
         od;
      od;
      Add (gens, newmat);
      Print(newmat, "\n");
    od;

    G := Source(rep);
    return GroupHomomorphismByImages(G, Group(gens), GeneratorsOfGroup(G), gens);
end;


decomvec := function(B, pat, v)
local u, k;
u := B*v;
k := 1;
#for c in pat do
#  M: = tomat(c, u{[k..k+c.dim*c.mul-1]});
#...
#od;
# ................   
end;
