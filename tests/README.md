# Rev3 unit tests

From the package root, run the three test entry points after adding YALMIP
and MOSEK to the MATLAB path and ensuring Optimization Toolbox is installed:

```matlab
run(fullfile('tests', 'run_forward_invariance_tests_rev3.m'))
run(fullfile('tests', 'run_sos_basis_tests_rev3.m'))
run(fullfile('tests', 'run_sos_postcheck_tests_rev3.m'))
```

`run_forward_invariance_tests_rev3.m` also needs `linprog` from the
Optimization Toolbox.  The regression test scans `T_sc = 0:0.1:1.5`; this
is intentionally broader than the paper configuration and checks only the
forward-invariance oracle.  It does not generate a paper certificate.

For the 17 quick tests used in the submission regression, exclude that
larger regression explicitly. From the package root:

```matlab
addpath(fullfile(pwd, 'data', 'double_integrator'))
addpath(fullfile(pwd, 'src', 'proposed', 'functions'))
addpath(fullfile(pwd, 'tests'))
yalmipRoot = fileparts(fileparts(which('sdpvar')));
addpath(fullfile(yalmipRoot, 'extras'), '-begin')
suite = matlab.unittest.TestSuite.fromClass(?ForwardInvarianceOracleRev3Test);
suite = suite(~contains({suite.Name}, 'realRev3SampleTimesRemainForwardInvariant'));
suite = [suite, matlab.unittest.TestSuite.fromClass(?SosBasisRev3Test), ...
    matlab.unittest.TestSuite.fromClass(?SosPostcheckGramBindingRev3Test)];
results = run(suite);
assert(all([results.Passed]), 'One or more quick tests failed.')
```
