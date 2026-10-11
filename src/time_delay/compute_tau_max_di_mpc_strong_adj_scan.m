%% Compute tau_max for double-integrator + explicit MPC (strong + adjacent scan)
% Method follows the workflow in example/compute_tau_max_theorem1_paper_sign_strong_adjacent_scan.m
% and uses a local generalized Theorem-1 feasibility checker for polytope regions.

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
releaseRoot = find_release_root(scriptDir);
addpath(scriptDir);
addpath(fullfile(releaseRoot, 'data', 'double_integrator'));
outputDir = ensure_generated_output_dir(releaseRoot);
configure_mosek_path();

solverName = string(getenv('TD_PWA_SOLVER'));
if strlength(strtrim(solverName)) == 0
    solverName = "sdpt3";
end
solverVerbose = 0;
mosekTol = str2double(getenv('TD_PWA_MOSEK_TOL'));
if ~isfinite(mosekTol) || mosekTol <= 0
    mosekTol = 1e-8;
end
mosekThreads = str2double(getenv('TD_PWA_MOSEK_THREADS'));
if ~isfinite(mosekThreads) || mosekThreads < 1
    mosekThreads = 1;
end
outputMatName = string(getenv('TD_PWA_OUTPUT_MAT'));
if strlength(strtrim(outputMatName)) == 0
    outputMatName = "double_integrator_tau_scan_strong_adjacent.mat";
end
summaryCsvName = string(getenv('TD_PWA_SUMMARY_CSV'));
if strlength(strtrim(summaryCsvName)) == 0
    summaryCsvName = "time_delay_paper_results.csv";
end
stateScale = parse_state_scale(getenv('TD_PWA_STATE_SCALE'));

if exist('theorem1_feasible_case_general_polytope', 'file') ~= 2
    error('Local generalized feasibility checker not found in: %s', scriptDir);
end

D = build_time_delay_data();

A = D.A;
a = D.a;
B = D.B;
K = D.K;
k = D.k;
E = D.E;
e = D.e;

M = numel(A);
nx = size(A{1}, 1);

% Normalize shapes to avoid row/column inconsistencies in LMI assembly.
B = B(:);
for i = 1:M
    K{i} = reshape(K{i}, 1, []);
    k{i} = k{i}(:);
    if numel(k{i}) ~= 1
        error('k{%d} must be scalar, got %d entries.', i, numel(k{i}));
    end
    k{i} = k{i}(1);

    a{i} = a{i}(:);
    e{i} = e{i}(:);
end

if ~isempty(stateScale)
    [A, a, B, K, E] = apply_state_scaling(A, a, B, K, E, stateScale);
    fprintf('Applied state scaling x = diag([%.6g %.6g]) z for LMI coordinates.\n', ...
        stateScale(1), stateScale(2));
end

fprintf('Loaded explicit-MPC closed-loop data: M=%d regions, nx=%d\n', M, nx);
fprintf('Plant conversion mode: %s\n', D.conversionMode);
fprintf('Input map B = [%.4f %.4f]^T\n', B(1), B(2));

%% Applicability checks (for this constructed case)
Aref = A{1};
maxADev = 0;
maxANorm = 0;
for i = 1:M
    maxADev = max(maxADev, norm(A{i} - Aref, inf));
    maxANorm = max(maxANorm, norm(a{i}, inf));
end
fprintf('Open-loop consistency check: max ||A_i - A_1||_inf = %.3e\n', maxADev);
fprintf('Affine-term check: max ||a_i||_inf = %.3e\n', maxANorm);

%% Auto-classify I0 via cancellation ||a_i + B*k_i||_inf <= tol
% Keep the same criterion used in the previous adaptation work.
tolI0 = 1e-9;
[I0_auto, cancel_norms] = classify_I0_from_cancellation(a, B, k, tolI0);
fprintf('I0 (auto): ');
fprintf('%d ', I0_auto);
fprintf('\n');
for i = 1:M
    fprintf('  region %d: ||a_i + B*k_i||_inf = %.3e\n', i, cancel_norms(i));
end

%% Solver/options and search settings
opts = make_yalmip_solver_options(solverName, solverVerbose, mosekTol, mosekThreads);
enforceLambdaEntrywise = true;

searchCfg.tauLowInit = 0.01;
searchCfg.tauHighInit = 0.10;
searchCfg.tauCap = 5.00;
searchCfg.tol = 1e-4;
searchCfg.maxIter = 50;
searchCfg.tauLowFloor = 1e-6;
searchCfg.thetaCandidateEps = 1e-3;

% Strong + adjacent scan settings (same style as reference script).
thetaUpperScanList = [0.99 0.97 0.95 0.93 0.90];
Nmargin = 1e-3;
strongGammaLowerBound = 1e-5;

%% Delta bounds
[deltaGlobal, deltaAdjacent, adjPairs] = compute_delta_bounds_general(K, k, E, e, nx, opts);

fprintf('Pairwise/global DeltaK = %.6f, Deltak = %.6f\n', deltaGlobal.DeltaK, deltaGlobal.Deltak);
fprintf('Adjacent DeltaK = %.6f, Deltak = %.6f, adjacent pairs = %d\n', ...
    deltaAdjacent.DeltaK, deltaAdjacent.Deltak, size(adjPairs,1));

fprintf('Adjacent pairs by closure-touch test:\n');
for n = 1:size(adjPairs,1)
    fprintf('  (%d, %d)\n', adjPairs(n,1), adjPairs(n,2));
end

%% Theta-target scan
scanResults = scan_theta_upper_targets(thetaUpperScanList, searchCfg.thetaCandidateEps, deltaAdjacent, Nmargin, ...
    A, a, B, K, k, E, e, I0_auto, nx, opts, enforceLambdaEntrywise, strongGammaLowerBound, searchCfg);

fprintf('\nStrong + adjacent theta scan results:\n');
for s = 1:numel(scanResults)
    r = scanResults(s);
    if r.invalidTheta
        fprintf('  thetaUpper=%.3f -> invalid thetaDesign=%.6f\n', r.thetaUpper, r.thetaDesign);
        continue;
    end

    if r.lowFeasibleFound
        lowFeasStr = sprintf('%.6f', r.lowFeasibleTau);
    else
        lowFeasStr = sprintf('none(<=%.1e)', r.lowSearchFloor);
    end

    fprintf(['  thetaUpper=%.3f, thetaDesign=%.6f, N=%.6f, tauBest=%.4f, ' ...
        'lowFeasible=%s, feasible=%d\n'], ...
        r.thetaUpper, r.thetaDesign, r.N_enforce, r.tauBest, lowFeasStr, r.lowFeasibleFound);
end

% Best over scanned theta targets.
validIdx = find(~[scanResults.invalidTheta]);
if isempty(validIdx)
    bestIdx = NaN;
    bestTau = NaN;
else
    tauVals = [scanResults(validIdx).tauBest];
    [bestTau, loc] = max(tauVals);
    bestIdx = validIdx(loc);
end

if ~isnan(bestIdx)
    fprintf('\nBest tau over scanned thetaUpper list: %.6f at thetaUpper=%.3f\n', ...
        bestTau, scanResults(bestIdx).thetaUpper);
end

results = struct();
results.sourceData = D.description;
results.I0_auto = I0_auto;
results.cancel_norms = cancel_norms;
results.maxADeviation = maxADev;
results.maxANorm = maxANorm;
results.deltaGlobal = deltaGlobal;
results.deltaAdjacent = deltaAdjacent;
results.adjPairs = adjPairs;
results.thetaUpperScanList = thetaUpperScanList;
results.scanResults = scanResults;
results.strongGammaLowerBound = strongGammaLowerBound;
results.Nmargin = Nmargin;
results.searchCfg = searchCfg;
results.solver = solverName;
results.mosekTol = mosekTol;
results.mosekThreads = mosekThreads;
results.stateScale = stateScale;
results.bestTau = bestTau;
results.bestScanIndex = bestIdx;

if strcmp(getenv('AMPC_SAVE_MAT'), '1')
    save(fullfile(outputDir, char(outputMatName)), 'results');
    fprintf('Saved optional MAT: %s\n', outputMatName);
end

if strlength(strtrim(summaryCsvName)) > 0
    summary = struct2table(scanResults);
    summary.solver = repmat(solverName, height(summary), 1);
    summary.mosekTol = repmat(mosekTol, height(summary), 1);
    summary.mosekThreads = repmat(mosekThreads, height(summary), 1);
    if isempty(stateScale)
        summary.stateScale_x1 = NaN(height(summary), 1);
        summary.stateScale_x2 = NaN(height(summary), 1);
    else
        summary.stateScale_x1 = repmat(stateScale(1), height(summary), 1);
        summary.stateScale_x2 = repmat(stateScale(2), height(summary), 1);
    end
    summary.bestTauOverall = repmat(bestTau, height(summary), 1);
    writetable(summary, fullfile(outputDir, char(summaryCsvName)));
    fprintf('Saved summary CSV: %s\n', summaryCsvName);
end


function scanResults = scan_theta_upper_targets(thetaUpperList, thetaEps, deltaAdjacent, Nmargin, ...
    A, a, B, K, k, E, e, I0, nx, opts, enforceLambdaEntrywise, gammaLowerBound, cfg)

n = numel(thetaUpperList);
scanResults = repmat(struct( ...
    'thetaUpper', NaN, 'thetaDesign', NaN, 'invalidTheta', false, ...
    'N_required', NaN, 'N_enforce', NaN, ...
    'tauBest', NaN, 'lowFeasibleTau', NaN, 'lowFeasibleFound', false, ...
    'lowInitTau', cfg.tauLowInit, 'lowSearchFloor', cfg.tauLowFloor, ...
    'problem', NaN), n, 1);

for s = 1:n
    thUp = thetaUpperList(s);
    thDesign = thUp - thetaEps;

    scanResults(s).thetaUpper = thUp;
    scanResults(s).thetaDesign = thDesign;

    if thDesign <= 0 || thDesign >= 1
        scanResults(s).invalidTheta = true;
        continue;
    end

    Nreq = (deltaAdjacent.DeltaK^2) / thDesign;
    Nenf = (1 + Nmargin) * Nreq;
    scanResults(s).N_required = Nreq;
    scanResults(s).N_enforce = Nenf;

    [tauBest, info, cert] = tau_max_bisection( ...
        A, a, B, K, k, E, e, I0, nx, opts, true, enforceLambdaEntrywise, Nenf, gammaLowerBound, cfg);

    scanResults(s).tauBest = tauBest;
    scanResults(s).lowFeasibleTau = info.lowFeasibleTau;
    scanResults(s).lowFeasibleFound = info.lowFeasibleFound;
    scanResults(s).lowInitTau = info.lowInitTau;
    scanResults(s).lowSearchFloor = info.lowSearchFloor;
    if isfield(cert, 'problem')
        scanResults(s).problem = cert.problem;
    end
end
end


function [tauBest, info, certBest] = tau_max_bisection(A, a, B, K, k, E, e, I0, nx, opts, ...
    enforceEtaGamma, enforceLambdaEntrywise, etaGammaScale, gammaLowerBound, cfg)

isFeasible = @(tau) theorem1_feasible_case_general_polytope( ...
    tau, A, a, B, K, k, E, e, I0, nx, opts, enforceEtaGamma, enforceLambdaEntrywise, etaGammaScale, gammaLowerBound);

lo = cfg.tauLowInit;
hi = cfg.tauHighInit;
lowFloor = cfg.tauLowFloor;
certBest = struct();

% Ensure lower bracket is feasible.
while lo > lowFloor && ~isFeasible(lo)
    lo = 0.5 * lo;
end
lowFeasibleFound = isFeasible(lo);
if lowFeasibleFound
    lowFeasibleTau = lo;
else
    lowFeasibleTau = NaN;
end

if ~lowFeasibleFound
    tauBest = -1;
    info = struct('lo', NaN, 'hi', NaN, 'iter', 0, 'isLowerBoundOnly', false, ...
        'lowFeasibleTau', lowFeasibleTau, 'lowFeasibleFound', false, ...
        'lowInitTau', cfg.tauLowInit, 'lowSearchFloor', lowFloor);
    certBest.problem = -1;
    certBest.info = 'No feasible point found near zero.';
    certBest.etaGammaScale = etaGammaScale;
    certBest.gammaLowerBound = gammaLowerBound;
    certBest.P = [];
    certBest.R = [];
    certBest.X = [];
    certBest.eta = [];
    certBest.gamma = [];
    return;
end

% Ensure upper bracket is infeasible, expanding if needed.
while hi < cfg.tauCap && isFeasible(hi)
    lo = hi;
    hi = min(1.5 * hi, cfg.tauCap);
    if hi == cfg.tauCap
        break;
    end
end

if hi == cfg.tauCap && isFeasible(hi)
    tauBest = hi;
    info = struct('lo', lo, 'hi', hi, 'iter', 0, 'isLowerBoundOnly', true, ...
        'lowFeasibleTau', lowFeasibleTau, 'lowFeasibleFound', true, ...
        'lowInitTau', cfg.tauLowInit, 'lowSearchFloor', lowFloor);
    [~, certBest] = theorem1_feasible_case_general_polytope( ...
        tauBest, A, a, B, K, k, E, e, I0, nx, opts, enforceEtaGamma, enforceLambdaEntrywise, etaGammaScale, gammaLowerBound);
    return;
end

iter = 0;
while (hi - lo) > cfg.tol && iter < cfg.maxIter
    mid = 0.5 * (lo + hi);
    if isFeasible(mid)
        lo = mid;
    else
        hi = mid;
    end
    iter = iter + 1;
end

tauBest = lo;
info = struct('lo', lo, 'hi', hi, 'iter', iter, 'isLowerBoundOnly', false, ...
    'lowFeasibleTau', lowFeasibleTau, 'lowFeasibleFound', true, ...
    'lowInitTau', cfg.tauLowInit, 'lowSearchFloor', lowFloor);
[~, certBest] = theorem1_feasible_case_general_polytope( ...
    tauBest, A, a, B, K, k, E, e, I0, nx, opts, enforceEtaGamma, enforceLambdaEntrywise, etaGammaScale, gammaLowerBound);
end


function [I0, norms] = classify_I0_from_cancellation(a, B, k, tol)
M = numel(a);
norms = zeros(M,1);
I0 = [];
for i = 1:M
    ci = a{i} + B*k{i};
    norms(i) = norm(ci, inf);
    if norms(i) <= tol
        I0(end+1) = i; %#ok<AGROW>
    end
end
end


function [deltaGlobal, deltaAdjacent, adjPairs] = compute_delta_bounds_general(K, k, E, e, nx, optsFeas)
M = numel(K);

DeltaK_global = 0;
Deltak_global = 0;
for i = 1:M
    for j = 1:M
        DeltaK_global = max(DeltaK_global, norm(K{i} - K{j}, 2));
        Deltak_global = max(Deltak_global, norm(k{i} - k{j}, 2));
    end
end

adjPairs = find_adjacent_pairs_by_closure_touch(E, e, nx, optsFeas);

DeltaK_adj = 0;
Deltak_adj = 0;
for n = 1:size(adjPairs,1)
    i = adjPairs(n,1);
    j = adjPairs(n,2);
    DeltaK_adj = max(DeltaK_adj, norm(K{i} - K{j}, 2));
    Deltak_adj = max(Deltak_adj, norm(k{i} - k{j}, 2));
end

% Fallback: if no adjacent pairs were found numerically, use global values.
if isempty(adjPairs)
    DeltaK_adj = DeltaK_global;
    Deltak_adj = Deltak_global;
end

deltaGlobal = struct('DeltaK', DeltaK_global, 'Deltak', Deltak_global);
deltaAdjacent = struct('DeltaK', DeltaK_adj, 'Deltak', Deltak_adj);
end


function adjPairs = find_adjacent_pairs_by_closure_touch(E, e, nx, optsFeas)
M = numel(E);
adjPairs = [];
tol = 1e-8;

for i = 1:M
    for j = i+1:M
        if closures_intersect(E{i}, e{i}, E{j}, e{j}, nx, tol, optsFeas)
            adjPairs = [adjPairs; i j]; %#ok<AGROW>
        end
    end
end
end


function stateScale = parse_state_scale(raw)
raw = string(raw);
if strlength(strtrim(raw)) == 0
    stateScale = [];
    return;
end
parts = regexp(char(raw), '[,;\\s]+', 'split');
parts = parts(~cellfun('isempty', parts));
vals = str2double(parts);
if numel(vals) ~= 2 || any(~isfinite(vals)) || any(vals <= 0)
    error('TD_PWA_STATE_SCALE must contain two positive numbers, e.g., "20,5".');
end
stateScale = vals(:);
end


function [A, a, B, K, E] = apply_state_scaling(A, a, B, K, E, stateScale)
Sx = diag(stateScale(:));
SxInv = diag(1 ./ stateScale(:));
for ii = 1:numel(A)
    A{ii} = SxInv * A{ii} * Sx;
    a{ii} = SxInv * a{ii};
    K{ii} = K{ii} * Sx;
    E{ii} = E{ii} * Sx;
end
B = SxInv * B;
end


function opts = make_yalmip_solver_options(solverName, verbose, mosekTol, mosekThreads)
opts = sdpsettings('solver', char(solverName), 'verbose', verbose);
opts.cachesolvers = 1;
if strcmpi(solverName, "mosek")
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = mosekTol;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = mosekTol;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = mosekTol;
    opts.mosek.MSK_IPAR_NUM_THREADS = mosekThreads;
end
end


function tf = closures_intersect(Ei, ei, Ej, ej, nx, tol, opts)
x = sdpvar(nx,1);
C = [Ei*x + ei >= -tol, Ej*x + ej >= -tol];
sol = optimize(C, [], opts);
tf = (sol.problem == 0);
end


function root = find_release_root(start_dir)
root = start_dir;
while true
    if exist(fullfile(root, 'data', 'double_integrator', 'answer_regions.m'), 'file') == 2
        return;
    end
    parent = fileparts(root);
    if strcmp(parent, root) || isempty(parent)
        error('Could not find release root from %s.', start_dir);
    end
    root = parent;
end
end


function output_dir = ensure_generated_output_dir(release_root)
output_dir = fullfile(release_root, 'results', 'generated');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
end


function configure_mosek_path()
mosek_path = getenv('AMPC_MOSEK_PATH');
if isempty(mosek_path) || ~isfolder(mosek_path)
    return;
end
if ~any(strcmpi(strsplit(path, pathsep), mosek_path))
    addpath(mosek_path);
end
end
