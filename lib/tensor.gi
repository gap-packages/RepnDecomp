InstallGlobalFunction( TensorProductOfMatricesObj, function(Fam, A, B)
    # check dimensions match the family dimensions
    if DimensionsMat(A) <> Fam!.dimA or DimensionsMat(B) <> Fam!.dimB then
        Error("Dimensions of tensor product don't match family!");
    fi;

    return Objectify(NewType(Fam, IsTensorProductOfMatricesObj and IsTensorProductPairRep),
                     [A, B]);
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
               IsIdenticalObj,
               [IsTensorProductOfMatricesObj and IsTensorProductPairRep,
                IsTensorProductOfMatricesObj and IsTensorProductPairRep],
               function(x, y)
                   return TensorProductOfMatricesObj(FamilyObj(x),
                                                     x![1]*y![1],
                                                     x![2]*y![2]);
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
