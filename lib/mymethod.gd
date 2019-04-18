#! @Chapter Computing decompositions of representations

#! @Section Original methods

#! @Arguments rho

#! @Returns A record

#! @Description Calculates some values given by other functions in
#! this package, but without summing over G (in case G is very
#! large). Most of the time, G is a symmetric group (or nearly one) so
#! calculation of conjugacy classes is easy.

#! <A>irreps</A> is the complete list of irreps involved in the direct
#! sum decomposition of <A>rho</A>, this can be given in case the
#! default (running Dixon's algorithm) is too expensive, or e.g. you
#! don't want representations over Cyclotomics.

#! If you have an orthonormal basis for the centraliser of <A>rho</A>,
#! you can pass it in as <A>rho_cent_basis</A> and it will be used to
#! speed up calculations.

#! The return value of this function is a record with fields:

#! * basis: same as BlockDiagonalBasisOfRepresentation

#! * diagonal_rep: same as BlockDiagonalRepresentation

#! * decomposition: same as IrreducibleDecompositionCollected

#! * centralizer_basis: same as RepresentationCentralizerBlocks

#! When I say "the same", I mean up to reordering and isomorphism as
#! representations.
DeclareAttribute( "REPN_ComputeUsingMyMethod", IsGroupHomomorphism );

#! @Arguments rho

#! @Returns A record in the same format as <Ref
#! Func="REPN_ComputeUsingMyMethod" />.

#! @Description Does the same thing as <Ref
#! Func="REPN_ComputeUsingMyMethod" />, but first splits the
#! representation into canonical summands. This might reduce the size
#! of the matrices we need to work over, so could much faster.
DeclareAttribute( "REPN_ComputeUsingMyMethodCanonical", IsGroupHomomorphism );
