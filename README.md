# NiBundleAdjustment

Generate Jacobians for Bundle Adjustment application with NiLang and ForwardDiff.

The motivation is to beat the benchmark in this [paper](https://arxiv.org/abs/1807.10129), for the glory of Julia community!

[![Build Status](https://travis-ci.com/JuliaReverse/NiBundleAdjustment.jl.svg?branch=master)](https://travis-ci.com/JuliaReverse/NiBundleAdjustment.jl)
[![Codecov](https://codecov.io/gh/JuliaReverse/NiBundleAdjustment.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaReverse/NiBundleAdjustment.jl)

## Get started!

Open a Julia REPL and type `]` to enter `pkg` mode and then type
```julia pkg
pkg> dev git@github.com:JuliaReverse/NiBundleAdjustment.jl.git
pkg> add ForwardDiff BenchmarkTools
```

Then in a bash shell, type the following commands to open the benchmark file in Atom.
```bash
$ julia ~/.julia/dev/NiBundleAdjustment/benchmarks/benchmark.jl
```

You will see results like:
```julia repl
Normal Objective
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     60.111 ns (0.00% GC)
  median time:      60.332 ns (0.00% GC)
  mean time:        61.376 ns (0.00% GC)
  maximum time:     162.091 ns (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     982
Reversible Objective
BenchmarkTools.Trial: 
  memory estimate:  144 bytes
  allocs estimate:  3
  --------------
  minimum time:     151.439 ns (0.00% GC)
  median time:      152.798 ns (0.00% GC)
  mean time:        171.360 ns (6.59% GC)
  maximum time:     5.764 μs (96.40% GC)
  --------------
  samples:          10000
  evals/sample:     818
NiLang Gradient
BenchmarkTools.Trial: 
  memory estimate:  332.28 MiB
  allocs estimate:  6543109
  --------------
  minimum time:     278.000 ms (7.06% GC)
  median time:      391.476 ms (24.87% GC)
  mean time:        394.913 ms (25.61% GC)
  maximum time:     488.938 ms (28.86% GC)
  --------------
  samples:          13
  evals/sample:     1ForwardDiff Gradient
BenchmarkTools.Trial: 
  memory estimate:  528.84 MiB
  allocs estimate:  4907333
  --------------
  minimum time:     591.759 ms (6.08% GC)
  median time:      901.409 ms (16.11% GC)
  mean time:        882.480 ms (17.03% GC)
  maximum time:     1.141 s (27.03% GC)
  --------------
  samples:          6
  evals/sample:     1
```

Note: the memory usage is not yet fully optimized. We can still trade some space for time.

It corresponds to the second column of ADBench paper
![ADBench](benchmarks/adbench.png)

We see our ForwardDiff result is 16x faster than the original version.
NiLang is even faster. Here the Jacobian is computed by glueing the Jacobians of a function with 15 input parameters and 2 output parameters. This is why NiLang does not have much advantage comparing with forward mode AD.
