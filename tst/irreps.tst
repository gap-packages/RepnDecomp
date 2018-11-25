gap> Read("tst/utils.g");
gap> G:=AlternatingGroup(5);;
gap> P:=PermutationGModule(G,Rationals);;
gap> h:=GroupHomomorphismByImages(G,Group(P.generators));;
gap> h:=ConvertRhoIfNeeded@RepnDecomp(h);;
gap> l:=DecomposeRepresentationIrreducible(h); # check the dimensions are correct
[ <vector space over Cyclotomics, with 1 generators>,
  <vector space over Cyclotomics, with 4 generators> ]
gap> ForAll(l, space -> IsGInvariant(h, space)); # check the spaces are actually subrepresentations
true
gap> G:=SymmetricGroup(3);;
gap> R:=RegularActionHomomorphism(G);;
gap> h:=GroupHomomorphismByImages(G,Image(R,G));;
gap> h:=ConvertRhoIfNeeded@RepnDecomp(h);;
gap> l:=DecomposeRepresentationIrreducible(h); # check the dimensions are correct
[ <vector space over Cyclotomics, with 1 generators>,
  <vector space over Cyclotomics, with 1 generators>,
  <vector space over Cyclotomics, with 2 generators>,
  <vector space over Cyclotomics, with 2 generators> ]
gap> ForAll(l, space -> IsGInvariant(h, space)); # check the spaces are actually subrepresentations
true
