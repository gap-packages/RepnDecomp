gap> BlockDiagonalMatrix([]); # no blocks
[  ]
gap> BlockDiagonalMatrix([[[1]], [], [[1,0],[0,1]]]); # some empty blocks
[ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ]
gap> # some differently sized nonempty blocks
gap> BlockDiagonalMatrix([IdentityMat(2), IdentityMat(3), IdentityMat(10)]) = IdentityMat(15);
true
