# model/ local notes

Loaded automatically when working on files in this directory. Full map:
`docs/agent/ARCHITECTURE.md`. File-by-file index: `docs/agent/FILE_INDEX.md`.

## Dead files in this directory, do not read or edit

- `modelnonrev.cpp` and `modelnonrev.h` are 0 bytes and absent from `CMakeLists.txt`.
- `modelgtr.cpp` is absent from `CMakeLists.txt` and includes `modelgtr.h`, which does not
  exist, so it cannot compile. Its contents are superseded by `modelmarkov.cpp`.

## Build

`CMakeLists.txt` in this directory is a hand-maintained explicit source list with no globbing.
A new `.cpp` or `.h` here that is not added to it is silently never compiled, and the symptom is
a link error rather than a missing-file error.

## Invariants that are easy to break

- Optimizer parameter vectors are 1-indexed, `variables[1..ndim]`, and `variables[0]` is unused.
  Rate arrays are 0-indexed. Every `memcpy` against a variable vector is `variables+1`.
- `targetFunk` returns the negative log-likelihood, because the optimizers minimize.
- `getVariables` must return whether anything actually changed. Callers skip recomputation when
  it returns false.
- After any change to Q or to `state_freq`, call `decomposeRateMatrix()` and then
  `phylo_tree->clearAllPartialLH()`. Omitting the second gives wrong likelihoods with no warning.
- Arrays the likelihood kernels read must be allocated with `aligned_alloc<double>` or
  `ensure_aligned_allocated`, never plain `new double[]`, or the vectorized kernels crash.
- `ModelMarkov::setRates()` is `ASSERT(0)`. Override it if the model uses an indirect
  parameterization, and actually call it from `getVariables`.
- Complex eigen arrays `ceval`, `cevec`, `cinv_evec` alias the same allocation as `eigenvalues`,
  `eigenvectors`, `inv_eigenvectors`. Do not free or reallocate one without the other.
- `PartitionModel::targetFunk` runs partitions under OpenMP, so anything it touches must be
  thread-safe or per-partition.

## Conventions

Use `outError`, `outWarning`, and `ASSERT` rather than `throw`, `assert`, or `std::cerr`. Gate
diagnostics on `verbose_mode >= VB_MED`. Read global options through `Params::getInstance()`
rather than threading new constructor arguments. Match the indentation of the file you are
editing, since older files use tabs and newer ones use four spaces.

## Reversible and non-reversible

One boolean, `is_reversible`, separates the two. For non-reversible models `rates[]` holds all
`n(n-1)` off-diagonal entries of Q in row-major order, the tree is rooted, the eigen spectrum is
complex, and `state_freq` is an output derived from Q on every evaluation rather than an input.
