# Some useful filters I need to define
InstallMethod( IsFiniteGroupLinearRepresentation,
               "checks if rho is a linear rep of finite group",
               [ IsGroupHomomorphism ],
               rho -> IsGroupHomomorphism(rho) and IsFinite(Source(rho)) and IsMatrixGroup(Range(rho)));

InstallMethod( IsFiniteGroupPermutationRepresentation,
               "finite to perm group",
               [ IsGroupHomomorphism ],
               rho -> IsGroupHomomorphism(rho) and IsFinite(Source(rho)) and IsPermGroup(Range(rho)));
