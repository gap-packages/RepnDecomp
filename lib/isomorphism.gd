#! @Chapter Isomorphisms between representations

#! @Section Finding explicit isomorphisms

#! @Description Let $\rho : G \to GL(V)$ and $\tau : G \to GL(W)$. If
#! there exists a linear map $A : V \to W$ such that for all $g \in
#! G$, $\tau(g)A = A\rho(g)$, this function returns one such $A$. $A$
#! is the isomorphism between the representations. If the
#! representations are not isomorphic, then fail is returned.
#!
#! There are three methods that we can use to compute an isomorphism
#! of linear representations, you can select one by passing options to
#! the function.
#!
#! * `use_kronecker`: Assumes the matrices are small enough that their
#!   Kronecker products can fit into memory. Uses <Ref
#!   Func="GroupSumBSGS" /> and `KroneckerProduct` to compute an
#!   element of the fixed subspace of $\rho \otimes \tau^*$.
#!
#! * `use_orbit_sum`: Finds an isomorphism by summing orbits of the
#!   the action of $\rho \otimes \tau^*$ on matrices. Note that orbits
#!   could be very large, so this could be as bad as summing over the
#!   whole group.
#!
#! * The default, sums over the whole group to compute the projection
#!   onto the fixed subspace.
#!
#! <P/>
#! @InsertChunk Example_LinearRepresentationIsomorphism
#! <P/>
#! @Arguments rho, tau[, rho_cent_basis, tau_cent_basis]
#! @Returns A matrix $A$ or fail
DeclareGlobalFunction( "LinearRepresentationIsomorphism" );


#! @Description Gives the same result as <Ref
#! Func="LinearRepresentationIsomorphism" />, but this function uses a
#! simpler method which always involves summing over $G$, without
#! using <Ref Func="GroupSumBSGS" />. This might be useful in some
#! cases if computing a good BSGS is difficult. However, for all cases
#! that have been tested, it is slow (as the name suggests).
#!
#! <P/>
#! @InsertChunk Example_LinearRepresentationIsomorphismSlow
#! <P/>
#! @Arguments rho, tau
#! @Returns A matrix $A$ or fail
DeclareGlobalFunction( "LinearRepresentationIsomorphismSlow" );


#! @Section Testing isomorphisms

#! @Description Since representations of finite groups over $\mathbb{C}$ are
#! determined by their characters, it is easy to check whether two
#! representations are isomorphic by checking if they have the same
#! character. We try to use characters wherever possible.

#! <P/>
#! @InsertChunk Example_AreRepsIsomorphic
#! <P/>
#! @Arguments rho, tau

#! @Returns true if <A>rho</A> and <A>tau</A> are isomorphic as
#! representations, false otherwise.
DeclareGlobalFunction( "AreRepsIsomorphic" );

#! @Description This function tests if, for all $g \in G$, $A \rho(g)
#! = \tau(g) A$. That is, true is returned iff $A$ is the intertwining
#! operator taking $\rho$ to $\tau$.
#! that:

#! <P/>
#! @InsertChunk Example_IsLinearRepresentationIsomorphism
#! <P/>

#! @Arguments A, rho, tau
#! @Returns true if <A>rho</A> and <A>tau</A> are isomorphic as as
#! representations with the isomorphism given by the linear map
#! <A>A</A>
DeclareGlobalFunction( "IsLinearRepresentationIsomorphism" );
