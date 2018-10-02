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
