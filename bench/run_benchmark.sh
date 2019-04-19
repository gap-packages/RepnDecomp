#!/bin/bash

# number of examples to benchmark (generated randomly)
bench_times=30

write_bench_file() {
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
}

# finding intertwining operators has 3 methods, based on how you
# compute the image of the trivial projection: orbit sums, using
# kronecker/BSGS or naive summing, this is the only variable in the
# alternate method
write_bench_file "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps)" "mymethod_naive"
write_bench_file "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps, use_kronecker)" "mymethod_kronecker"
write_bench_file "REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps, use_orbit_sum)" "mymethod_orbit_sum"

# do the same for mymethodcanonical
write_bench_file "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps)" "mymethod_canonical_naive"
write_bench_file "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_kronecker)" "mymethod_canonical_kronecker"
write_bench_file "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_orbit_sum)" "mymethod_canonical_orbit_sum"

# and the same for mymethodcanonical in parallel, to see the
# difference
write_bench_file "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, parallel)" "mymethod_parallel_naive"
write_bench_file "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_kronecker, parallel)" "mymethod_parallel_kronecker"
write_bench_file "REPN_ComputeUsingMyMethodCanonical(rep.rep : irreps := rep.irreps, use_orbit_sum, parallel)" "mymethod_parallel_orbit_sum"

# serre's method has using kroneckers and not as the only variables we
# can really test randomly. need to get unitary reps and orthonormal
# centraliser bases to test cent_basis stuff, can't really do this
# randomly
write_bench_file "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps)" "serre_naive"
write_bench_file "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, use_kronecker)" "serre_kronecker"

# the same but in parallel
write_bench_file "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, parallel)" "serre_parallel_naive"
write_bench_file "REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps, use_kronecker, parallel)" "serre_parallel_kronecker"

# Runs the benchmarks for various functions, prints results
for bench in *.g; do
    echo "Running '$bench'"
    gap -q <$bench
done
