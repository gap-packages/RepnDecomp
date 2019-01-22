LoadPackage("RepnDecomp");;

# benchmarks a function to_bench that does something to
# representations, for a range of sizes of group and degree of
# representation
BenchRepForSmallGroups := function(to_bench, out)
    local size, num, id, G, irreps, rho, starttime, r, endtime;
    PrintTo(out);
    for size in [100, 110 .. 500] do
        # just testing ~10 groups should always be enough
        num := Minimum(10, NrSmallGroups(size));
        for id in [1..num] do
            G := SmallGroup(size, id);
            irreps := IrreducibleRepresentations(G);

            # the trivial rep
            rho := GroupHomomorphismByImages(G,
                                             Group(IdentityMat(1)),
                                             GeneratorsOfGroup(G),
                                             List(GeneratorsOfGroup(G), _ -> IdentityMat(1)));
            starttime := Runtime();
            r := to_bench(rho, irreps);
            endtime := Runtime();

            # group size, number of irreps/conjugacy classes, time taken
            AppendTo(out, Size(G), " ", Length(irreps), " ", endtime-starttime, "\n");
        od;
    od;
end;;

# benches function that does something to canonical summands
BenchCanonicalSummandForSmallGroups := function(to_bench, out)
    local size, num, id, G, irreps, rho, starttime, r, endtime, V;
    PrintTo(out);
    for size in [100, 110 .. 500] do
        # just testing ~10 groups should always be enough
        num := Minimum(10, NrSmallGroups(size));
        for id in [1..num] do
            G := SmallGroup(size, id);
            irreps := IrreducibleRepresentations(G);

            # so there's a fairly large canonical summand
            rho := DirectSumRepList(Replicate@RepnDecomp(irreps[1], 5));
            V := Cyclotomics^5;
            starttime := Runtime();
            r := to_bench(rho, irreps[1], V);
            endtime := Runtime();

            # group size, number of irreps/conjugacy classes, time taken
            AppendTo(out, Size(G), " ", Length(irreps), " ", endtime-starttime, "\n");
        od;
    od;
end;;
