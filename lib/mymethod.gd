#! @Chapter Computing decompositions of representations

#! @Section Algorithms due to the authors

#! @Arguments rho
#!
#! @Returns A record in the same format as <Ref
#! Attr="REPN_ComputeUsingSerre" Label="for IsGroupHomomorphism" />
#!
#! @Description Computes the same values as <Ref
#! Attr="REPN_ComputeUsingSerre" Label="for IsGroupHomomorphism" />,
#! taking the same options. The heavy lifting of this method is done
#! by <Ref Func="LinearRepresentationIsomorphism" />, where there are
#! some further options that can be passed to influence algorithms
#! used.
DeclareAttribute( "REPN_ComputeUsingMyMethod", IsGroupHomomorphism );

#! @Arguments rho
#!
#! @Returns A record in the same format as <Ref
#! Attr="REPN_ComputeUsingMyMethod" Label="for IsGroupHomomorphism"
#! />.
#!
#! @Description Performs the same computation as <Ref
#! Attr="REPN_ComputeUsingMyMethod" Label="for IsGroupHomomorphism"
#! />, but first splits the representation into canonical summands
#! using <Ref Func="CanonicalDecomposition" />. This might reduce the
#! size of the matrices we need to work with significantly, so could
#! be much faster.
#!
#! If the option `parallel` is given, the decomposition of canonical
#! summands into irreps is done in parallel, which could be much
#! faster.
DeclareAttribute( "REPN_ComputeUsingMyMethodCanonical", IsGroupHomomorphism );
