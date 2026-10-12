# AMPC stability-verification code

This package accompanies the paper on a posteriori stability verification for
asynchronous MPC (AMPC). It contains the selected fixed-control-sampling-time
configurations used in the numerical study: the proposed Rev3 piecewise-
polynomial Lyapunov/SOS method and three adapted comparison approaches.

## Associated manuscript and citation

**SOS-Based Stability Verification for Asynchronous Execution of Constrained LQ-MPC**

Feng Xu, Shengyi Liu, Shubo Yang, and Xiangyang Wang

Institute for Aero Engine, Tsinghua University, Beijing 100084, China

As of 12 October 2026, the manuscript is prepared for submission to
*Automatica*. It has not yet been submitted to the journal, and no public
manuscript link or journal DOI is available. The paper link and citation
will be updated when a public manuscript link or publication details become
available.

Until then, the manuscript can be identified as follows:

```bibtex
@unpublished{Xu2026AMPCStability,
  author = {Xu, Feng and Liu, Shengyi and Yang, Shubo and Wang, Xiangyang},
  title = {{SOS}-Based Stability Verification for Asynchronous Execution of Constrained {LQ-MPC}},
  year = {2026},
  note = {Manuscript prepared for submission}
}
```

Code repository: [Asynchronous-LQ-MPC-Stability-Verification](https://github.com/2505733155/Asynchronous-LQ-MPC-Stability-Verification).

The `v1.0-submission` tag preserves the initial code snapshot of 11 October
2026. The `v1.0-github-preparation` tag adds this manuscript metadata on
12 October 2026; executable code and reference CSVs are identical to the
initial snapshot. These tags describe preparation snapshots, rather than
confirming a journal submission or acceptance.

## Scope and dependencies

The package requires MATLAB, YALMIP, and MOSEK. Its reference outputs were
generated and checked in the following software environment.

| Component | Validated version | Role |
|---|---|---|
| MATLAB | 24.2.0.3212159 (R2024b) Update 9 | Computing platform |
| YALMIP | R20250626 (26 June 2025) | SOS/LMI modeling |
| MOSEK MATLAB Toolbox and solver | 11.2.1 (build 2026-06-01) | SDP solver |
| Optimization Toolbox | 24.2 | Required by the proposed, PWQ, and KKT-lifted routines (`linprog` or `quadprog`) |
| Control System Toolbox | 24.2 | Optional: supplies `dare` for the KKT routine; the package otherwise uses its fixed-point fallback |
| Parallel Computing Toolbox | 24.2 | Optional: accelerates the Rev3 scan; the package otherwise runs serially |

The YALMIP version is the value returned by `yalmip('version')`. The validated
YALMIP release archive did not retain Git metadata, so R20250626 is recorded
as the precise available release identifier rather than an unavailable commit
hash. The package does not require the Model Predictive Control Toolbox.

Before running the paper cases, verify the installed stack in MATLAB:

```matlab
version
yalmip('version')
mosekopt('version')
```

A clean-path regression on 11 October 2026 also passed on MATLAB
24.2.0.3361896 (R2024b) Update 10, with the same YALMIP and MOSEK versions.
The public data validator, 17 quick unit tests, and all four paper entry
points completed successfully. The conclusion fields matched the committed
reference CSVs (10 proposed, 180 PWQ, five time-delay, and four KKT cases).
The larger forward-invariance-only regression was not rerun in this check.
The Update 9 environment above remains the provenance of the reference CSVs;
the Update 10 rerun does not replace them.

Other versions have not been validated against the committed reference CSVs;
SOS/SDP diagnostics and borderline solver outcomes can differ across numerical
software versions.

Install YALMIP and MOSEK separately; neither is redistributed here. MOSEK
should normally already be on the MATLAB path. If needed, set one optional
environment variable before starting MATLAB. Add the installed YALMIP release
to the MATLAB path as well:

```matlab
addpath(genpath('path-to-YALMIP'))
setenv('AMPC_MOSEK_PATH', 'path-to-the-MOSEK-MATLAB-toolbox')
```

From the package root, first validate the public controller data and then run
all four paper configurations:

```matlab
addpath(fullfile(pwd, 'data', 'double_integrator'))
validate_double_integrator_data
run_all_paper_cases
```

Individual entry points are in `src/proposed/`, `src/pwq/`,
`src/time_delay/`, and `src/kkt/`. Rev3 unit tests are documented in
[`tests/README.md`](tests/README.md).

## Configurations and reference results

| Configuration | Public entry point | Reference CSV |
|---|---|---|
| Proposed Rev3 piecewise-polynomial Lyapunov/SOS | `src/proposed/run_rev3_unbounded_final_physdense.m` | `results/reference/proposed_rev3_paper_results.csv` |
| PWQ Lyapunov SDP margin sweep | `src/pwq/run_pwq_rev1_mosek_paper_margin_sweep.m` | `results/reference/pwq_paper_margin_sweep_results.csv` |
| Sampled-data PWA time-delay LMI | `src/time_delay/run_time_delay_mosek_scaled_diag.m` | `results/reference/time_delay_paper_scaled20_5_results.csv` |
| KKT-lifted SOS global direct attempt | `src/kkt/run_kkt_lifted_sos_standard_mpc_global_scaled20_5_attempt.m` | `results/reference/kkt_global_paper_scaled20_5_results.csv` |

The proposed configuration uses `d_V=4`, `d_m=2`,
`epsilon=alpha=1e-2`, `T_sc=0.1:0.1:1.0`, a common state scaling
`diag(20,5)`, and MOSEK feasibility/gap tolerances of `1e-8`. The three
comparison entries use the selected paper settings stated in their wrappers;
all use the same `diag(20,5)` state scaling where the formulation contains
state coordinates.

The PWQ entry performs the complete paper sweep: 10 sampling times,
`epsilon=alpha=1e-2,1e-3,...,1e-10`, and the existing `scalar` and `matrix`
multiplier modes, for 180 feasibility attempts. The `1e-2` through `1e-9`
runs return no numerical candidate. At `1e-10`, the matrix-multiplier runs
return numerical candidates at all ten sampling times, whereas the scalar
multiplier runs retain their actual non-candidate outcome. Thus, the complete
CSV makes explicit the multiplier-mode qualification behind the paper's
statement about the additional `1e-10` run.

To reproduce either comparison independently from the package root, run:

```matlab
run(fullfile('src', 'pwq', 'run_pwq_rev1_mosek_paper_margin_sweep.m'))
run(fullfile('src', 'time_delay', 'run_time_delay_mosek_scaled_diag.m'))
```

The PWQ command writes the 180-row
`results/generated/pwq_paper_margin_sweep_results.csv`; the time-delay
command writes the five-row `results/generated/time_delay_paper_scaled20_5_results.csv`.
The complete PWQ sweep takes several minutes on the validated MATLAB/MOSEK
environment; `run_all_paper_cases` additionally runs the proposed and KKT
configurations.

Committed CSV files are reference outputs only. Reruns write CSVs and any
optional local artifacts to `results/generated/`, which is ignored by Git.
Set `AMPC_SAVE_MAT=1` only when a local MAT snapshot is wanted; MAT files are
not part of this release.

The time-delay entry scans `thetaUpper = [0.99, 0.97, 0.95, 0.93, 0.90]`,
with `thetaDesign = thetaUpper - 1e-3`. For each value, the lower feasibility
search begins at `0.01 s` and repeatedly halves the trial sampling-interval
upper bound until it reaches the configured `1e-6 s` floor. The five reference
rows all report `lowFeasibleFound=false`, `tauBest=-1`, and `problem=-1`;
this is the directly reproducible basis for the paper's statement that no
feasible numerical solution was returned even near an upper bound of
`1e-6 s`. Here `tauBest=-1` denotes that the script found no feasible lower
bracket, rather than a solver-certified instability conclusion.

The time-delay reference CSV was retained from the paper run. The public
package consolidates the controller at its higher stored precision, whereas
that archived run used rounded controller coefficients; therefore its
intermediate `N_required` and `N_enforce` values can differ by about 0.2% on a
current MATLAB/MOSEK installation. This data-precision effect does not change
the conclusion fields above.

The Update 10 regression reproduced the proposed residual/Gram diagnostics
and KKT diagnostics in the reference CSVs; their differing fields were run
times only. PWQ differences were confined to run times and non-conclusion
partition/transition diagnostics (LP witness coordinates, selected region
constraints, complement margin, and Chebyshev radii). The time-delay
`N_required` and `N_enforce` differences were approximately 0.17%, within the
data-precision effect stated above. All reference CSVs are retained unchanged.

## Random seeds

The KKT routine fixes the finite-point closed-loop diagnostic seed to `7`
(`cfg.collision_random_seed`). Its auxiliary random residual probes use
`rng(11 + meta.gram_dim + meta.degree)` for each block. These probes are
numerical diagnostics, not exact certificate validation. No explicit random
sampling was found in the other three paper entry points; their configured
point grids are deterministic. Fixed seeds do not eliminate solver/version
or parallel-runtime differences.

## Interpretation

The results concern fixed values of the control sampling time. They are
floating-point numerical candidates and associated diagnostics, not a
certificate over a continuous interval of sampling times or an exact rational
certificate. Solver infeasibility, an inconclusive diagnostic, or the lack of
a certificate for a baseline does **not** establish closed-loop instability.

The comparison matrix discussed in the paper also mentions a bounded-domain
PWA Lyapunov approach for approximate explicit MPC. It is not included here:
its bounded-analysis-domain scope differs from the selected study and no final
runnable implementation was used as a paper baseline. The package likewise
does not include internal TDO/PGD diagnostics.

## Layout

```text
data/double_integrator/   shared plant and seven-region PWA controller data
src/                      one implementation directory per included approach
tests/                    unit tests for the proposed route
results/reference/        four committed paper-reference CSV files
results/generated/        local rerun output, ignored by Git
```

This repository distributes only code and newly written documentation under
the MIT license. MATLAB, YALMIP, MOSEK, and cited literature remain subject
to their respective licenses.
