gap> # some random group
gap> G := SmallGroup(20, 4);;
gap> irreps := IrreducibleRepresentations(G);;
gap> # note deg irreps[2] = 1, deg irreps[5] = 2
gap> gens := GeneratorsOfGroup(G);;
gap> imgs := List(gens, g -> BlockDiagonalMatrix([Image(irreps[2], g), Image(irreps[2], g), Image(irreps[2], g), Image(irreps[5], g), Image(irreps[5], g)]));;
gap> rho := GroupHomomorphismByImages(G, Group(imgs), gens, imgs);;
gap> # rho should split into 3 lines and 2 planes
gap> IrreducibleDecompositionCollected(rho).decomp;
[ [  ],
  [
      rec( basis := [ [ 1, 0, 0, 0, 0, 0, 0 ] ],
          space := <vector space over Cyclotomics, with 1 generators> ),
      rec( basis := [ [ 0, 1, 0, 0, 0, 0, 0 ] ],
          space := <vector space over Cyclotomics, with 1 generators> ),
      rec( basis := [ [ 0, 0, 1, 0, 0, 0, 0 ] ],
          space := <vector space over Cyclotomics, with 1 generators> ) ], [  ],
  [  ],
  [ rec( basis := [ [ 0, 0, 0, 1, 0, 0, 0 ], [ 0, 0, 0, 0, 1, 0, 0 ] ],
          space := <vector space over Cyclotomics, with 2 generators> ),
      rec( basis := [ [ 0, 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 0, 1 ] ],
          space := <vector space over Cyclotomics, with 2 generators> ) ], [  ],
  [  ], [  ] ]