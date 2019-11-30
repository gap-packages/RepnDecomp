gap> V := VectorSpace(Rationals, [[1,0],[0,1]], [0,0]);;
gap> p := [[1,0],[0,0]];; # projects to a line
gap> W := MatrixImage@RepnDecomp(p, V);
<vector space over Rationals, with 2 generators>
gap> List(Basis(W));
[ [ 1, 0 ] ]
gap> Dimension(W);
1
gap> V := VectorSpace(Cyclotomics, [], 0);;
gap> MatrixImage@RepnDecomp(p, V); # do nothing if basis is empty
<vector space of dimension 0 over Cyclotomics>
