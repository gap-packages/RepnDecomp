gap> l := [1,2,3];;
gap> Take@RepnDecomp(l, 2);
[ 1, 2 ]
gap> Take@RepnDecomp(l, 3);
[ 1, 2, 3 ]
gap> Take@RepnDecomp(l, 0);
[  ]
gap> l; # check original is unchanged
[ 1, 2, 3 ]