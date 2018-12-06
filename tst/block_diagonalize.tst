gap> G:=SymmetricGroup(3);;
gap> R:=RegularActionHomomorphism(G);;
gap> h:=GroupHomomorphismByImages(G,Image(R,G));;
gap> h:=ConvertRhoIfNeeded@RepnDecomp(h);;
gap> f:=BlockDiagonalRepresentation(h);
[ (1,2,3), (1,2) ] ->
[ [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, E(3), 0, 0, 0 ],
      [ 0, 0, 0, E(3)^2, 0, 0 ], [ 0, 0, 0, 0, E(3), 0 ],
      [ 0, 0, 0, 0, 0, E(3)^2 ] ],
  [ [ 1, 0, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0, 0 ], [ 0, 0, 0, E(3)^2, 0, 0 ],
      [ 0, 0, E(3), 0, 0, 0 ], [ 0, 0, 0, 0, 0, E(3)^2 ],
      [ 0, 0, 0, 0, E(3), 0 ] ] ]