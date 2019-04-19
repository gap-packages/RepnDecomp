#! @Chapter Miscellaneous useful functions

#! @Section Efficient summing over groups

#! @Arguments G, summand

#! @Returns $\sum_{g \in G} \mbox{summand}(g)$

#! @Description Uses a basic stabiliser chain for $G$ to compute the
#! sum described above. This trick requires <A>summand</A> to be a
#! function (in the GAP sense) that defines a monoid homomorphism (in
#! the mathematical sense). The computation of the stabiliser chain
#! assumes <A>G</A> is a group.

#! More precisely, if we have the basic stabiliser chain:

#! $$\{1\} = G_1 \leq \ldots \leq G_n = G$$

#! We traverse the chain from $G_1$ to $G_n$, using the previous sum
#! $G_{i-1}$ to build the sum $G_i$. We do this by using the fact
#! that (writing $f$ for <A>summand</A>)

#! $$\sum_{g \in G_i} f(g) = \sum_{r_j} \left(\sum_{h \in G_{i-1}} f(h)\right) f(r_j)$$

#! where the $r_j$ are right coset representatives of $G_{i-1}$ in
#! $G_i$.

#! The condition on <A>summand</A> is satisfied if, for example, it is
#! a linear representation of a group <A>G</A>.
DeclareGlobalFunction( "GroupSumBSGS" );
