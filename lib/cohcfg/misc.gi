# Replace and extend GRAPE's OrbitNumbers function, using native GAP functions
# G is a group
# Omega is a set upon which G acts (please make sure it is closed under action by G, as this will not be checked!)
OrbitsInfo@ := function (G, Omega)
    local orbits, mapping, representatives,
        i, j;
    if IsPosInt(Omega) then Omega := [1..Omega]; # backwards-compatibility with GRAPE's OrbitNumbers() syntax
    elif not IsHomogeneousList(Omega) then
        Error("OrbitsInfo@: malformed Omega");
    fi;
    
    orbits := OrbitsDomain(G, Omega); # this returns a list of the orbits; WARNING: undefined behavior when Omega is not closed!
    mapping := EmptyPlist(Length(Omega));
    for i in [1..Length(orbits)] do
        for j in orbits[i] do
            mapping[j] := i; # element j in Omega is in the i-th orbit; WARNING: mapping domain may be too big when Omega is not closed!
        od;
    od;
    MakeImmutable(mapping);

    representatives := List(orbits, x -> x[1]); # retrieves the first element of each orbit
    MakeImmutable(representatives);
    
    return rec(
        orbits := orbits, # a list of the orbits of G acting on Omega
        mapping := mapping, # a mapping from Omega to orbits, i.e. mapping[i] = j iff i in Omega is in orbit j
            orbitNumbers := ~.mapping, # backwards-compatibility with GRAPE's OrbitNumbers() syntax
        representatives := representatives # a list whose i-th entry is a representative element from the i-th orbit of Omega
    );
end;
