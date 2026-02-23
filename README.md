# PMFRG for the XYZ model

## Examples

Refer to the `examples` directory and documentation therein 
for a minimal running example. 
See also the tests.


## Notes on the temporary state of the code

These notes might change as the code changes. 

### Versions of the code: T-flow and Lambda-flow

- T-flow and Lambda-flow are both implemented in  `src/PMFRG_xyz.jl`.

The goal of the current refactoring effort is to have both versions in the same codebase.

The differences between the two versions reside mainly 
in the functions `iS_`, `iG_`, and `iSKat_`.

### Testing

To run tests on the Temperature-flow version, which is the main code,
one can use the standard command
```
julia -t auto --project=. -e "using Pkg; Pkg.test()"
```

To run tests for the Lambda-flow version (a single end-to-end regression test):
```
julia -t auto --project=./test/regression/lambda_flow_regression/ ./test/regression/lambda_flow_regression/run_test.jl
```

The test suites still need to be merged.
