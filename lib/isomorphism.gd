#! @Chapter Isomorphisms between representations

#! @Section Finding explicit isomorphisms

#! @Arguments rho, tau

#! @Returns A matrix $A$ or fail

#! @Description Let $\rho : G \to GL(V)$ and $\tau : G \to GL(W)$. If
#! there exists a linear map $A : V \to W$ such that for all $g \in
#! G$, $\tau(g)A = A\rho(g)$, this function returns one such $A$. $A$
#! is the isomorphism between the representations. If the
#! representations are not isomorphic, then fail is returned.
DeclareGlobalFunction( "RepresentationIsomorphism" );

#! @Arguments rho, tau

#! @Returns true if <A>rho</A> and <A>tau</A> are isomorphic as
#! representations, false otherwise.

#! @Description Note that two representations are isomorphic iff they
#! give similar matrices, also iff they have the same irreducible
#! decomposition. We use characters to determine the latter: the first
#! is expensive to check for large degree representations.
DeclareGlobalFunction( "AreRepsIsomorphic" );
