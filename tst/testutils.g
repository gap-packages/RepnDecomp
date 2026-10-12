#
# Shared helpers for RepnDecomp's hand-written tests.
#

BindGlobal( "REPN_TEST_ConjugateRepresentation", function( rho, A )
    return ComposeHomFunction( rho, M -> A^-1 * M * A );
end );

BindGlobal( "REPN_TEST_RestrictedRepresentation", function( rho, summand )
    local basis;

    if IsRecord( summand ) and IsBound( summand.basis ) then
        basis := summand.basis;
    elif IsVectorSpace( summand ) then
        basis := List( Basis( summand ) );
    else
        Error( "test summand must be a record with a basis or a vector space" );
    fi;

    basis := Basis( VectorSpace( Cyclotomics, basis ) );
    return RestrictRep@RepnDecomp( rho, basis );
end );

BindGlobal( "REPN_TEST_IsIrreducible", function( rho )
    local chi;
    chi := CharacterOfRepresentation@RepnDecomp( rho );
    return InnerProductOfCharacters@RepnDecomp( chi, chi, Source( rho ) ) = 1;
end );

BindGlobal( "REPN_TEST_IsInvariantDecomposition", function( rho, decomposition )
    local bases, degree;

    bases := List( decomposition, function( summand )
        if IsRecord( summand ) and IsBound( summand.basis ) then
            return summand.basis;
        elif IsVectorSpace( summand ) then
            return List( Basis( summand ) );
        fi;
        Error( "test summand must be a record with a basis or a vector space" );
    end );
    degree := DegreeOfRepresentation( rho );

    return Sum( bases, Length ) = degree
           and RankMat( Concatenation( bases ) ) = degree
           and ForAll( bases,
                       basis -> IsGInvariant@RepnDecomp( rho, basis ) )
           and ForAll( decomposition,
                       summand -> REPN_TEST_IsIrreducible(
                           REPN_TEST_RestrictedRepresentation( rho, summand ) ) );
end );

BindGlobal( "REPN_TEST_IsCollectedDecomposition", function( rho, collected )
    local flat, restricted, i, j;

    flat := Flat( collected );
    if not REPN_TEST_IsInvariantDecomposition( rho, flat ) then
        return false;
    fi;

    restricted := List( collected,
                        family -> List( family,
                            summand -> REPN_TEST_RestrictedRepresentation(
                                rho, summand ) ) );

    for i in [ 1 .. Length( restricted ) ] do
        if not ForAll( restricted[i],
                       rep -> AreRepsIsomorphic( restricted[i][1], rep ) ) then
            return false;
        fi;

        if i < Length( restricted ) then
            for j in [ i + 1 .. Length( restricted ) ] do
                if AreRepsIsomorphic( restricted[i][1], restricted[j][1] ) then
                    return false;
                fi;
            od;
        fi;
    od;

    return true;
end );

BindGlobal( "REPN_TEST_CommutesWithRepresentation", function( rho, A )
    return ForAll( GeneratorsOfGroup( Source( rho ) ),
                   g -> A * Image( rho, g ) = Image( rho, g ) * A );
end );

BindGlobal( "REPN_TEST_IsBlockDiagonalization", function( rho, info )
    local P, degree;

    degree := DegreeOfRepresentation( rho );
    P := TransposedMat( info.basis );

    return Length( info.basis ) = degree
           and RankMat( info.basis ) = degree
           and ForAll( GeneratorsOfGroup( Source( rho ) ),
                       g -> Image( info.diagonal_rep, g )
                            = P^-1 * Image( rho, g ) * P );
end );

BindGlobal( "REPN_TEST_UpperUnitriangularMat", function( n )
    local A, i;

    A := IdentityMat( n );
    for i in [ 1 .. n - 1 ] do
        A[i][i + 1] := i;
    od;
    if n > 2 then
        A[1][n] := -1;
    fi;
    return A;
end );
