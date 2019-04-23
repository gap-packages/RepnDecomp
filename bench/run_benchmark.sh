#!/bin/bash

# first we test some known examples, like:

# * the defining representations of some symmetric groups
# * the regular representations of same
# * some graph automorphism groups (Johnson, Hamming)

# we can test really big symmetric groups since the defining
# representation is quite small
run_symmetric_test () {
    bench_str="$1"
    name="$2"
    cat <<EOF > ${name}.g
LoadPackage("RepnDecomp");;
rhos := List([1..10], n -> ConvertRhoIfNeeded@RepnDecomp(IdentityMapping(SymmetricGroup(n))));
BenchList@RepnDecomp(rep -> ${bench_str}, rhos, "${name}.csv");
EOF
    echo Running ${name}.g
    gap -q < ${name}.g
}

# can only test for small symmetric groups, regular repn is really big
run_regular_test () {
    bench_str="$1"
    name="$2"
    cat <<EOF > ${name}.g
LoadPackage("RepnDecomp");;
rhos := List([1..6], n -> ConvertRhoIfNeeded@RepnDecomp(RegularActionHomomorphism(SymmetricGroup(n))));
BenchList@RepnDecomp(rep -> ${bench_str}, rhos, "${name}.csv");
EOF
    echo Running ${name}.g
    gap -q < ${name}.g
}

# TODO: actually do the tests, decide which to test etc

# next, we do the random tests

# number of examples to benchmark (generated randomly)
bench_times=30

run_random_test() {
    bench_str="$1"
    name="$2"
    cat <<EOF > ${name}.g
LoadPackage("RepnDecomp");;
# options for the random generator
opt := rec(lo := 50,
           hi := 200,
           num_irreps := 2,
           min_multiplicity := 1,
           max_multiplicity := 2,
           max_total_degree := 10,
           restrict_small_degree := false,
           small_degree := 10);;
BenchMany@RepnDecomp(rep -> ${bench_str}, "${name}.csv", ${bench_times}, opt);
EOF
    echo Running ${name}.g
    gap -q < ${name}.g
}

# finding intertwining operators has 3 methods, based on how you
# compute the image of the trivial projection: orbit sums, using
# kronecker/BSGS or naive summing, this is the only variable in the
# alternate method
run_random_test "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps)" "mymethod_naive"
run_random_test "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps, use_kronecker)" "mymethod_kronecker"
run_random_test "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps, use_orbit_sum)" "mymethod_orbit_sum"

# do the same for mymethodcanonical
run_random_test "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps)" "mymethod_canonical_naive"
run_random_test "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_kronecker)" "mymethod_canonical_kronecker"
run_random_test "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_orbit_sum)" "mymethod_canonical_orbit_sum"

# and the same for mymethodcanonical in parallel, to see the
# difference
run_random_test "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, parallel)" "mymethod_parallel_naive"
run_random_test "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_kronecker, parallel)" "mymethod_parallel_kronecker"
run_random_test "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_orbit_sum, parallel)" "mymethod_parallel_orbit_sum"

# serre's method has using kroneckers and not as the only variables we
# can really test randomly. need to get unitary reps and orthonormal
# centraliser bases to test cent_basis stuff, can't really do this
# randomly
run_random_test "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps)" "serre_naive"
run_random_test "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, use_kronecker)" "serre_kronecker"

# the same but in parallel
run_random_test "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, parallel)" "serre_parallel_naive"
run_random_test "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, use_kronecker, parallel)" "serre_parallel_kronecker"
