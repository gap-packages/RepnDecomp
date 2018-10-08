gap> G:=SymmetricGroup(3);
Sym( [ 1 .. 3 ] )
gap> R:=RegularActionHomomorphism(G);
<action isomorphism>
gap> h:=GroupHomomorphismByImages(G,Image(R,G));
[ (1,2,3), (1,2) ] -> [ (1,4,5)(2,3,6), (1,3)(2,4)(5,6) ]
gap> f:=BlockDiagonalizeRepresentation(h);
[ (1,2,3), (1,2) ] -> 
[ [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, E(3), 0, 0, 0 ], 
      [ 0, 0, 0, E(3)^2, 0, 0 ], [ 0, 0, 0, 0, E(3), 0 ], 
      [ 0, 0, 0, 0, 0, E(3)^2 ] ], 
  [ [ 1, 0, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0, 0 ], [ 0, 0, 0, E(3)^2, 0, 0 ], 
      [ 0, 0, E(3), 0, 0, 0 ], [ 0, 0, 0, 0, 0, E(3)^2 ], 
      [ 0, 0, 0, 0, E(3), 0 ] ] ]
gap> 