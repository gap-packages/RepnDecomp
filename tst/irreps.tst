gap> G:=AlternatingGroup(5);
Alt( [ 1 .. 5 ] )
gap> P:=PermutationGModule(G,Rationals);
rec( IsOverFiniteField := false, dimension := 5, field := Rationals, 
  generators := 
    [ 
      [ [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ], 
          [ 0, 0, 0, 0, 1 ], [ 1, 0, 0, 0, 0 ] ], 
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 0, 1, 0 ], 
          [ 0, 0, 0, 0, 1 ], [ 0, 0, 1, 0, 0 ] ] ], isMTXModule := true )
gap> h:=GroupHomomorphismByImages(G,Group(P.generators));
[ (1,2,3,4,5), (3,4,5) ] -> 
[ 
  [ [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ],
      [ 1, 0, 0, 0, 0 ] ], 
  [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 0, 1, 0 ], 
      [ 0, 0, 0, 0, 1 ], [ 0, 0, 1, 0, 0 ] ] ]
gap> l:=DecomposeRepresentationIrreducible(h);
[ <vector space over Cyclotomics, with 1 generators>, 
  <vector space over Cyclotomics, with 4 generators> ]
gap> G:=SymmetricGroup(3);
Sym( [ 1 .. 3 ] )
gap> R:=RegularActionHomomorphism(G);
<action isomorphism>
gap> h:=GroupHomomorphismByImages(G,Image(R,G));
[ (1,2,3), (1,2) ] -> [ (1,4,5)(2,3,6), (1,3)(2,4)(5,6) ]
gap> DecomposeRepresentationIrreducible(h);
[ <vector space over Cyclotomics, with 1 generators>, 
  <vector space over Cyclotomics, with 1 generators>, 
  <vector space over Cyclotomics, with 2 generators>, 
  <vector space over Cyclotomics, with 2 generators> ]
