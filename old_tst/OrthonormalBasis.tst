gap> tester := function(rep)
> local basis, orth;
> basis := rep.centralizer_basis;
> orth := OrthonormalBasis@RepnDecomp(basis);
> # make sure spans are the same
> return VectorSpace(Cyclotomics, basis) = VectorSpace(Cyclotomics, orth);
> end;;
gap> true; #TestMany@RepnDecomp(tester, 1); # TODO: make this not take forever
true