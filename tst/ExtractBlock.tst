gap> mat := IdentityMat(4);
[ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ]
gap> ExtractBlock@RepnDecomp(mat, 1, 1, 2);
[ [ 1, 0 ], [ 0, 1 ] ]
gap> ExtractBlock@RepnDecomp(mat, 1, 2, 2);
[ [ 0, 0 ], [ 0, 0 ] ]
gap> ExtractBlock@RepnDecomp(mat, 2, 2, 2);
[ [ 1, 0 ], [ 0, 1 ] ]
gap> ExtractBlock@RepnDecomp(mat, 1, 1, 3);
[ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ]