gap> l := [1,2,3];;
gap> L := Replicate@RepnDecomp(l, 3);
[ [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ]
gap> L[1][1] := 2;;
gap> L; # check only L[1] changes
[ [ 2, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ]
gap> l; # check original is unchanged
[ 1, 2, 3 ]
gap> Replicate@RepnDecomp(l, 0);
[  ]
gap> Replicate@RepnDecomp(3, 3);
[ 3, 3, 3 ]