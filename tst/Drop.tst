gap> l := [1,2,3];;
gap> Drop@RepnDecomp(l, 1);
[ 2, 3 ]
gap> Drop@RepnDecomp(l, 0);
[ 1, 2, 3 ]
gap> Drop@RepnDecomp(l, 3);
[  ]
gap> l; # check original is unchanged
[ 1, 2, 3 ]