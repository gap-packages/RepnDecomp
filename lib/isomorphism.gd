#! @Chapter Isomorphisms between representations

#! @Section Finding explicit isomorphisms

#! @Arguments rho, tau[, rho_cent_basis, tau_cent_basis]

#! @Returns A matrix $A$ or fail

#! @Description Let $\rho : G \to GL(V)$ and $\tau : G \to GL(W)$. If
#! there exists a linear map $A : V \to W$ such that for all $g \in
#! G$, $\tau(g)A = A\rho(g)$, this function returns one such $A$. $A$
#! is the isomorphism between the representations. If the
#! representations are not isomorphic, then fail is returned.

#! If no bases for the centralizers are given, we have no option but
#! to sum over $G$ to perform this calculation.  If bases for the
#! centralizers are given, we can use them to calculate class/group
#! sums and thus avoid summing over $G$.
DeclareGlobalFunction( "LinearRepresentationIsomorphism" );

#! @Arguments rho, tau

#! @Returns A matrix $A$ or fail

#! @Description The same as <Ref
#! Func="LinearRepresentationIsomorphism" />, but this function uses a
#! simpler method which always involves summing over $G$. This is slow
#! for large $G$, but might be fast in the special case of a very
#! large group and very small degree representation.
DeclareGlobalFunction( "LinearRepresentationIsomorphismSlow" );

#! @Arguments rho, tau

#! @Returns true if <A>rho</A> and <A>tau</A> are isomorphic as
#! representations, false otherwise.

#! @Description Note that two representations are isomorphic iff they
#! give similar matrices, also iff they have the same irreducible
#! decomposition. We use characters to determine the latter: the first
#! is expensive to check for large degree representations.
DeclareGlobalFunction( "AreRepsIsomorphic" );
