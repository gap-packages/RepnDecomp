gap> groups := List([2..4], SymmetricGroup);;
gap> replists := List(groups, IrreducibleRepresentations);;
gap> # check we get the right number (product) of irreps
gap> Length(TensorProductRepLists(replists[1], replists[2])) = Length(IrreducibleRepresentations(DirectProduct(groups[1], groups[2])));
true
gap> Length(TensorProductRepLists(replists[2], replists[3])) = Length(IrreducibleRepresentations(DirectProduct(groups[2], groups[3])));
true
gap> rep1 := GroupHomomorphismByFunction(groups[1], Group(IdentityMat(2)), g -> IdentityMat(2));;
gap> rep2 := GroupHomomorphismByFunction(groups[3], Group(IdentityMat(3)), g -> IdentityMat(3));;
gap> prod := TensorProductRepLists([rep1], [rep2]);
[ [ (1,2), (3,4,5,6), (3,4) ] ->
    [ [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ],
          [ 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 1 ] ],
      [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ],
          [ 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 1 ] ],
      [ [ 1, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ],
          [ 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 1 ] ] ] ]
gap> Source(prod[1]) = DirectProduct(groups[1], groups[3]);
true
gap> Image(prod[1], (1,2)) = IdentityMat(6);
true