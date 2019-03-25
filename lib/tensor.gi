InstallGlobalFunction( TensorProductOfMatricesObj, function(Fam, A, B)
    # check dimensions match the family dimensions
    if DimensionsMat(A) <> Fam!.dimA or DimensionsMat(B) <> Fam!.dimB then
        Error("Dimensions of tensor product don't match family!");
    fi;

    return Objectify(NewType(Fam, IsTensorProductOfMatricesObj and IsTensorProductPairRep),
                     [A, B]);
end );

InstallGlobalFunction( TensorProductOfMatrices, function(A, B)
    local F;
    F := TensorProductOfMatricesFamily(DimensionsMat(A), DimensionsMat(B));
    return TensorProductOfMatricesObj(F, A, B);
end );

IsDimensions@ := function(x)
    return IsList(x) and Length(x) = 2 and ForAll(x, IsPosInt);
end;

InstallGlobalFunction( TensorProductOfMatricesFamily, function(dimA, dimB)
    local F, G;

    if not IsDimensions@(dimA) or not IsDimensions@(dimB) then
        Error("<dimA> and <dimB> must be pairs of positive integers");
    fi;

    F := NewFamily( Concatenation("TensorProductOfMatrices", String(dimA), "x", String(dimB)),
                    IsTensorProductOfMatricesObj );

    F!.dimA := dimA;
    F!.dimB := dimB;

    return F;
end );

InstallMethod( \*,
               "for tensor products of matrices",
               [IsTensorProductOfMatricesObj and IsTensorProductPairRep,
                IsTensorProductOfMatricesObj and IsTensorProductPairRep],
               function(x, y)
                   return TensorProductOfMatricesObj(FamilyObj(x),
                                                     x![1]*y![1],
                                                     x![2]*y![2]);
               end );

InstallMethod( \*,
               "for applying tensor products to matrices",
               [IsTensorProductOfMatricesObj and IsTensorProductPairRep,
                IsMatrix],
               function(prod, A)
                   # A pure tensor product x tensor y acts on e_i
                   # tensor e_j, giving x(e_i) tensor y(e_j). The
                   # action extends to all matrices by linearity.
                   #
                   # Linearity also lets us avoid materialising the
                   # Kronecker product, but we have to spend O(n^4)
                   # time looping over indices
                   local i, j, result, e, n, m, ith_col, jth_col;

                   n := Length(A);
                   m := Length(A[1]);

                   # check the sizes are all correct
                   if n <> Length(prod![1][1]) or m <> Length(prod![2][1]) then
                       Error("Tensor product is wrong size to act on <A>!");
                   fi;

                   result := NullMat(n, m);

                   for i in [1..n] do
                       for j in [1..m] do
                           ith_col := TransposedMat([TransposedMat(prod![1])[i]]);
                           jth_col := [TransposedMat(prod![2])[j]];
                           result := result + A[i][j] * KroneckerProduct(ith_col, jth_col);
                       od;
                   od;

                   return result;
               end );

InstallMethod( \=,
               "for tensor products of matrices",
               [IsTensorProductOfMatricesObj and IsTensorProductPairRep,
                IsTensorProductOfMatricesObj and IsTensorProductPairRep],
               function(x, y)
                   return x![1] = y![1] and x![2] = y![2];
               end );

InstallMethod( OneOp,
               "for tensor products of matrices",
               [IsTensorProductOfMatricesObj and IsTensorProductPairRep],
               function(x)
                   local F;
                   F := FamilyObj(x);
                   return TensorProductOfMatricesObj(F,
                                                     IdentityMat(F!.dimA[1]),
                                                     IdentityMat(F!.dimB[1]));
               end );

InstallMethod( InverseOp,
               "for tensor products of matrices",
               [IsTensorProductOfMatricesObj and IsTensorProductPairRep],
               function(x)
                   return TensorProductOfMatricesObj(FamilyObj(x),
                                                     x![1]^-1,
                                                     x![2]^-1);
               end );

# We don't want to have to print the kronecker product matrix, would
# use too much memory.
InstallMethod( PrintObj,
               "for tensor products of matrices",
               [IsTensorProductOfMatricesObj and IsTensorProductPairRep],
               function(x)
                   Print(x![1], " tensor ", x![2]);
               end );

InstallMethod( IsGeneratorsOfMagmaWithInverses,
               "for tensor products of matrices",
               [CategoryCollections(IsTensorProductOfMatricesObj)],
               function(coll)
                   # just check if all matrices are invertible
                   return ForAll(coll, function(tensor)
                                    if IsTensorProductPairRep(tensor) then
                                        return IsGeneratorsOfMagmaWithInverses([tensor![1]])
                                               and IsGeneratorsOfMagmaWithInverses([tensor![2]]);
                                    elif IsTensorProductKroneckerRep(tensor) then
                                        return IsGeneratorsOfMagmaWithInverses([tensor![1]]);
                                    else
                                        return fail;
                                    fi;
                                end );
               end );

# takes tensor product of two matrix representations
InstallGlobalFunction( TensorProductOfRepresentations, function(rho, tau)
    local F, G, deg_rho, deg_tau;

    deg_rho := DegreeOfRepresentation(rho);
    deg_tau := DegreeOfRepresentation(tau);
    F := TensorProductOfMatricesFamily([deg_rho, deg_rho], [deg_tau, deg_tau]);

    if Source(rho) <> Source(tau) then
        Error("Not representations of the same group!");
    fi;

    G := Source(rho);

    return FuncToHom@(G, g -> TensorProductOfMatricesObj(F, Image(rho, g), Image(tau, g)));
end );

# this is the old way, sometimes useful if you need to sum
InstallGlobalFunction( KroneckerProductOfRepresentations, function(rho, tau)
    return FuncToHom@(Source(rho), g -> KroneckerProduct(Image(rho, g), Image(tau, g)));
end );

# don't have to materialise matrices if you just want the character
InstallGlobalFunction( CharacterOfTensorProductOfRepresentations, function(rho)
    return function(g)
        local im;
        im := Image(rho, g);
        return Trace(im![1]) * Trace(im![2]);
    end;
end );
