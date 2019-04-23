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
rhos := List([2..6], n -> ConvertRhoIfNeeded@RepnDecomp(IdentityMapping(SymmetricGroup(n))));;
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
rhos := List([2..6], n -> ConvertRhoIfNeeded@RepnDecomp(RegularActionHomomorphism(SymmetricGroup(n))));;
BenchList@RepnDecomp(rep -> ${bench_str}, rhos, "${name}.csv");
EOF
    echo Running ${name}.g
    gap -q < ${name}.g
}

# next, we do the random tests

# number of examples to benchmark (generated randomly)
bench_times=30

run_random_test() {
    bench_str="$1"
    name="$2"
    cat <<EOF > ${name}.g
LoadPackage("RepnDecomp");;
# options for the random generator
opt := rec(lo := 10,
           hi := 100,
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

# runs scripts for all combinations we want to test, using the given generation script
run_all_tests() {
    fn=$1
    tag=$2
    # finding intertwining operators has 3 methods, based on how you
    # compute the image of the trivial projection: orbit sums, using
    # kronecker/BSGS or naive summing, this is the only variable in the
    # alternate method
    $fn "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps)" "${tag}_mymethod_naive"
    $fn "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps, use_kronecker)" "${tag}_mymethod_kronecker"
    $fn "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps, use_orbit_sum)" "${tag}_mymethod_orbit_sum"

    # do the same for mymethodcanonical
    $fn "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps)" "${tag}_mymethod_canonical_naive"
    $fn "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_kronecker)" "${tag}_mymethod_canonical_kronecker"
    $fn "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_orbit_sum)" "${tag}_mymethod_canonical_orbit_sum"

    # and the same for mymethodcanonical in parallel, to see the
    # difference
    $fn "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, parallel)" "${tag}_mymethod_parallel_naive"
    $fn "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_kronecker, parallel)" "${tag}_mymethod_parallel_kronecker"
    $fn "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_orbit_sum, parallel)" "${tag}_mymethod_parallel_orbit_sum"

    # serre's method has using kroneckers and not as the only variables we
    # can really test randomly. need to get unitary reps and orthonormal
    # centraliser bases to test cent_basis stuff, can't really do this
    # randomly
    $fn "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps)" "${tag}_serre_naive"
    $fn "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, use_kronecker)" "${tag}_serre_kronecker"

    # the same but in parallel
    $fn "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, parallel)" "${tag}_serre_parallel_naive"
    $fn "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, use_kronecker, parallel)" "${tag}_serre_parallel_kronecker"
}

run_all_tests run_symmetric_test "symmetric"
run_all_tests run_regular_test "regular"
#run_all_tests run_random_test "random"
