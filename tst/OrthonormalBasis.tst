gap> tester := function(rep)
> local basis, orth;
> basis := rep.centralizer_basis;
> orth := OrthonormalBasis@RepnDecomp(basis);
> # make sure spans are the same
> return VectorSpace(Cyclotomics, basis) = VectorSpace(Cyclotomics, orth);
> end;;
gap> TestMany@RepnDecomp(tester, 2);
true