#ReadPackage("cohcfg", "lib/classmatr.g");
# this file is currently undocumented

GetSubCC@ := function (CC, H) # CC a CohCfg object, H a subgroup
    local SubCC, reps,
        i, x;

    # Sanity check
    if not IsCohCfg@(CC) then
        Error("GetSubCC@: CC must be an IsCohCfg@ object.\n");
    elif not IsSubgroup(CC.group, H) then
        Error("GetSubCC@: H must be a subgroup of the underlying group of CC.\n");
    fi;
    
    # Create sub-coherent configuration
    SubCC := CohCfgFromPermGroup@(H, CC.Omega);
    SubCC.isSubCC := true;

    # Determine which G-orbital each H-orbital is contained in
    reps := CCTransversal@(SubCC);
    SubCC.parentOrbitalMapping := List([1..SubCC.dimension], x -> CCElementContainingPair@(CC, reps[x]));
    SubCC.childOrbitalMapping := List([1..CC.dimension], x -> Positions(SubCC.parentOrbitalMapping, x));
    SubCC.parentCC := CC;

    return SubCC;
    # SubCC contains all members of a CohCfg object record, as well as the following:
    #   .isSubCC - identifies the record as a SubCohCfg object
    #   .parentCC - a handle for the parent CohCfg object from which SubCC was created
    #   .parentOrbitalMapping - a list mapping SubCC's own orbitals to those of the parent CC
    #   .childOrbitalMapping - a list mapping orbitals of the parent CC to lists of orbitals of SubCC
end;

IsSubCC@ := function (CC) # CC an object
    return IsCohCfg@ and IsBound(CC.isSubCC) and CC.isSubCC;
end;

GetParentCC@ := function (CC) # CC a SubCC object
    # Sanity check
    if not IsSubCC@(CC) then
        Error("GetParentCC@: CC must be a sub- coherent configuration object.\n");
    fi;

    # Get parent CC object
    return CC.parentCC;
end;
