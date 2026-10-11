%% Rev3 fixed-Ts explicit-MPC PWA certificate via piecewise polynomial SOS
% Formal certificate generator for:
%   Piecewise polynomial Lyapunov + SOS / higher-order S-procedure
% Here "higher-order" means higher-degree SOS multipliers, not generator products.
%
% Default run scans:
%   Ts = 0:0.1:1.5, d_V = 4, d_mult = 2
% using YALMIP/MOSEK.

clear; clc;

cfg = make_cfg();
script_dir = get_script_dir();
release_root = find_release_root(script_dir);
addpath(script_dir);
addpath(fullfile(script_dir, 'functions'));
addpath(fullfile(release_root, 'data', 'double_integrator'));
ensure_yalmip_export_precedence();
configure_mosek_path(cfg);
ensure_matlab_linprog_precedence();

output_dir = ensure_generated_output_dir(release_root);
out_csv = fullfile(output_dir, char(cfg.output_stem + "_results.csv"));
out_mat = fullfile(output_dir, char(cfg.output_stem + "_certificate.mat"));

fprintf('Rev3 piecewise-polynomial SOS certificate generator\n');
fprintf('Output CSV: %s\n', out_csv);
fprintf('Optional MAT output: %s\n', out_mat);
fprintf('Run label: %s\n', cfg.run_label);

reg = answer_regions();
scale_ctx = build_scale_context(reg, cfg);

rows = [];
cert_bundle = struct();
cases = build_case_list(cfg);

if should_use_parallel(cfg, numel(cases))
    pool = start_parallel_pool(cfg);
    fprintf('Parallel Ts scan enabled with %d worker(s).\n', pool.NumWorkers);
    row_cells = cell(numel(cases), 1);
    cert_cells = cell(numel(cases), 1);
    parfor kk = 1:numel(cases)
        setup_parallel_worker(release_root, script_dir, cfg);
        this_case = cases(kk);
        fprintf('\n[Rev3/parfor] Ts=%.12g, d_V=%d, d_mult=%d, eps=%.3g, alpha=%.3g\n', ...
            this_case.Ts, this_case.dV, this_case.dmult, cfg.epsilon, cfg.alpha);
        [row_cells{kk}, cert_cells{kk}] = solve_one_case( ...
            this_case.Ts, this_case.dV, this_case.dmult, reg, scale_ctx, cfg);
    end

    for kk = 1:numel(cases)
        rows = append_row(rows, row_cells{kk});
        cert_bundle.(cases(kk).key) = cert_cells{kk};
        print_case_summary(row_cells{kk});
    end
else
    fprintf('Parallel Ts scan disabled; running serial scan.\n');
    for kk = 1:numel(cases)
        this_case = cases(kk);
        fprintf('\n[Rev3] Ts=%.12g, d_V=%d, d_mult=%d, eps=%.3g, alpha=%.3g\n', ...
            this_case.Ts, this_case.dV, this_case.dmult, cfg.epsilon, cfg.alpha);
        [row, cert] = solve_one_case(this_case.Ts, this_case.dV, this_case.dmult, reg, scale_ctx, cfg);
        rows = append_row(rows, row);
        cert_bundle.(this_case.key) = cert;
        print_case_summary(row);
    end
end

T = struct2table(rows);
writetable(T, out_csv);
if strcmp(getenv('AMPC_SAVE_MAT'), '1')
    save(out_mat, 'cert_bundle', 'cfg', 'scale_ctx');
end

disp(T(:, {'Ts','d_V','d_mult','certificate_status','solver_problem_code', ...
    'reachable_transition_count','sos_constraint_count','postcheck_min_primal_residual','solve_time_sec'}));


function cfg = make_cfg()
cfg.run_label = "scan_Ts_0_0p1_1p5_dV4_dmult2_eps1em6_alpha1em6";
cfg.output_stem = "region_piecewise_polynomial_sos_rev3";
cfg.sample_times = 0:0.1:1.5;
cfg.dV_values = 4;
cfg.dmult_values = 2;
cfg.use_parallel = true;
cfg.parallel_workers = 4;

% To run the main Rev3 attempt, change to:
%   cfg.run_label = "main_scan_dV4_dmult2";
%   cfg.sample_times = 0:0.1:1.5;
%   cfg.dmult_values = 2;

cfg.epsilon = 1e-6;
cfg.alpha = 1e-6;
cfg.state_scale = [20, 5];
cfg.solver = "mosek";
cfg.mosek_path = getenv("AMPC_MOSEK_PATH");
cfg.verbose = 0;
cfg.yalmip_sos_model = 0;
cfg.yalmip_sos_scale = 1;
% Explicit candidate monomials are passed to solvesos, so YALMIP Newton
% monomial reduction is bypassed by construction.
cfg.yalmip_sos_newton = 0;
cfg.mosek_primal_feas_tol = 1e-8;
cfg.mosek_dual_feas_tol = 1e-8;
cfg.mosek_rel_gap_tol = 1e-8;
cfg.mosek_max_iter = 200;
cfg.mosek_num_threads = 1;
cfg.lp_tol = 1e-9;
cfg.forward_invariance_tol = 1e-6;
cfg.postcheck_tol = 1e-6;
cfg.sanity_tol = 1e-6;
cfg.gram_eig_tol = 1e-7;
cfg.coeff_match_tol = 1e-2;
cfg.sos_oracle_tol_factor = 10;
cfg.zero_row_tol = 1e-11;
cfg.redundancy_tol = 1e-9;
cfg.near_degenerate_radius_tol = 1e-8;
cfg.partition_tol = 1e-9;
cfg.origin_equilibrium_tol = 1e-10;
cfg.origin_structure_tol = 1e-8;

cfg.normalization_mode = "fixed_margin_plus_box";
cfg.coefficient_bound = 1e3;
cfg.max_sos_constraints = 500;
cfg.max_total_gram_dim = 2500;
cfg.grid_x1 = linspace(-30, 30, 61);
cfg.grid_x2 = linspace(-12, 12, 49);
cfg.dense_postcheck_enabled = true;
cfg.dense_postcheck_max_Ts = 1.0;
cfg.dense_grid_coordinate_system = "scaled_state";
cfg.dense_grid_scaled_x1 = linspace(-30, 30, 81);
cfg.dense_grid_scaled_x2 = linspace(-12, 12, 65);

cfg = apply_environment_overrides(cfg);
end


function cfg = apply_environment_overrides(cfg)
cfg.run_label = override_string(cfg.run_label, "REV3_RUN_LABEL");
cfg.output_stem = override_string(cfg.output_stem, "REV3_OUTPUT_STEM");
cfg.epsilon = override_number(cfg.epsilon, "REV3_EPSILON");
cfg.alpha = override_number(cfg.alpha, "REV3_ALPHA");
cfg.mosek_primal_feas_tol = override_number(cfg.mosek_primal_feas_tol, "REV3_MOSEK_PRIMAL_FEAS_TOL");
cfg.mosek_dual_feas_tol = override_number(cfg.mosek_dual_feas_tol, "REV3_MOSEK_DUAL_FEAS_TOL");
cfg.mosek_rel_gap_tol = override_number(cfg.mosek_rel_gap_tol, "REV3_MOSEK_REL_GAP_TOL");
shared_mosek_tol = getenv("REV3_MOSEK_TOL");
if strlength(string(shared_mosek_tol)) > 0
    tol = str2double(shared_mosek_tol);
    if ~isfinite(tol) || tol <= 0
        error('Invalid REV3_MOSEK_TOL: %s', shared_mosek_tol);
    end
    cfg.mosek_primal_feas_tol = tol;
    cfg.mosek_dual_feas_tol = tol;
    cfg.mosek_rel_gap_tol = tol;
end
cfg.use_parallel = override_logical(cfg.use_parallel, "REV3_USE_PARALLEL");
cfg.parallel_workers = override_number(cfg.parallel_workers, "REV3_PARALLEL_WORKERS");
cfg.sample_times = override_numeric_vector(cfg.sample_times, "REV3_SAMPLE_TIMES");
cfg = override_dense_physical_box(cfg);
end


function value = override_string(value, env_name)
env_value = getenv(char(env_name));
if strlength(string(env_value)) > 0
    value = string(env_value);
end
end


function value = override_number(value, env_name)
env_value = getenv(char(env_name));
if strlength(string(env_value)) == 0
    return;
end
new_value = str2double(env_value);
if ~isfinite(new_value)
    error('Invalid %s: %s', char(env_name), env_value);
end
value = new_value;
end


function value = override_logical(value, env_name)
env_value = lower(strtrim(string(getenv(char(env_name)))));
if strlength(env_value) == 0
    return;
end
if any(env_value == ["1", "true", "yes", "on"])
    value = true;
elseif any(env_value == ["0", "false", "no", "off"])
    value = false;
else
    error('Invalid %s: %s', char(env_name), env_value);
end
end


function value = override_numeric_vector(value, env_name)
env_value = strtrim(string(getenv(char(env_name))));
if strlength(env_value) == 0
    return;
end
tokens = regexp(char(env_value), '^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$', 'tokens', 'once');
if ~isempty(tokens)
    a = str2double(tokens{1});
    step = str2double(tokens{2});
    b = str2double(tokens{3});
    value = a:step:b;
    return;
end
parts = regexp(char(env_value), '[,;\s]+', 'split');
parts = parts(~cellfun('isempty', parts));
new_value = str2double(parts);
if isempty(new_value) || any(~isfinite(new_value))
    error('Invalid %s: %s', char(env_name), env_value);
end
value = new_value;
end


function cfg = override_dense_physical_box(cfg)
box_value = strtrim(string(getenv("REV3_DENSE_PHYSICAL_BOX")));
if strlength(box_value) == 0
    return;
end

box_parts = regexp(char(box_value), '[,;\s]+', 'split');
box_parts = box_parts(~cellfun('isempty', box_parts));
box = str2double(box_parts);
if numel(box) ~= 4 || any(~isfinite(box))
    error('REV3_DENSE_PHYSICAL_BOX must contain four finite numbers: x1_min x1_max x2_min x2_max.');
end
if box(1) >= box(2) || box(3) >= box(4)
    error('REV3_DENSE_PHYSICAL_BOX must have min < max in each coordinate.');
end

count_value = strtrim(string(getenv("REV3_DENSE_GRID_COUNTS")));
if strlength(count_value) == 0
    counts = [numel(cfg.dense_grid_scaled_x1), numel(cfg.dense_grid_scaled_x2)];
else
    count_parts = regexp(char(count_value), '[,;\s]+', 'split');
    count_parts = count_parts(~cellfun('isempty', count_parts));
    counts = round(str2double(count_parts));
    if numel(counts) ~= 2 || any(~isfinite(counts)) || any(counts < 2)
        error('REV3_DENSE_GRID_COUNTS must contain two integers >= 2.');
    end
end

cfg.dense_grid_coordinate_system = "scaled_state_from_physical_box";
cfg.dense_grid_scaled_x1 = linspace(box(1) / cfg.state_scale(1), box(2) / cfg.state_scale(1), counts(1));
cfg.dense_grid_scaled_x2 = linspace(box(3) / cfg.state_scale(2), box(4) / cfg.state_scale(2), counts(2));
end


function cases = build_case_list(cfg)
cases = struct('Ts', {}, 'dV', {}, 'dmult', {}, 'key', {});
kk = 0;
for tt = 1:numel(cfg.sample_times)
    for dd = 1:numel(cfg.dV_values)
        for mm = 1:numel(cfg.dmult_values)
            kk = kk + 1;
            Ts = cfg.sample_times(tt);
            dV = cfg.dV_values(dd);
            dmult = cfg.dmult_values(mm);
            cases(kk).Ts = Ts; %#ok<AGROW>
            cases(kk).dV = dV; %#ok<AGROW>
            cases(kk).dmult = dmult; %#ok<AGROW>
            cases(kk).key = matlab.lang.makeValidName(sprintf('Ts_%s_dV_%d_dmult_%d', ...
                num2str(Ts), dV, dmult)); %#ok<AGROW>
        end
    end
end
end


function tf = should_use_parallel(cfg, case_count)
tf = false;
if case_count <= 1 || ~isfield(cfg, 'use_parallel') || ~cfg.use_parallel
    return;
end
if isempty(which('parpool')) || isempty(which('gcp'))
    warning('Parallel Computing Toolbox not available. Falling back to serial scan.');
    return;
end
tf = true;
end


function pool = start_parallel_pool(cfg)
pool = gcp('nocreate');
if isempty(pool)
    workers = max(1, round(cfg.parallel_workers));
    pool = parpool('local', workers);
elseif pool.NumWorkers ~= cfg.parallel_workers
    warning('Using existing parallel pool with %d worker(s), cfg.parallel_workers=%d.', ...
        pool.NumWorkers, cfg.parallel_workers);
end
end


function setup_parallel_worker(release_root, script_dir, cfg)
addpath(script_dir);
addpath(fullfile(script_dir, 'functions'));
addpath(fullfile(release_root, 'data', 'double_integrator'));
ensure_yalmip_export_precedence();
configure_mosek_path(cfg);
ensure_matlab_linprog_precedence();
end


function print_case_summary(row)
fprintf('  Ts=%.12g status=%s, transitions=%d, sos_constraints=%d, min_check=%.4g, solve_time=%.2fs\n', ...
    row.Ts, row.certificate_status, row.reachable_transition_count, ...
    row.sos_constraint_count, row.postcheck_min_primal_residual, ...
    row.solve_time_sec);
end


function [row, cert] = solve_one_case(Ts, dV, dmult, reg_law, scale_ctx, cfg)
t_case = tic;
row = base_row(Ts, dV, dmult, cfg);
cert = empty_cert(Ts, dV, dmult, cfg);

if isempty(which('sdpvar')) || isempty(which('sos')) || isempty(which('solvesos')) || ...
        isempty(which('polynomial')) || isempty(which('coefficients'))
    row.certificate_status = "solver_failed";
    row.solver_info = "yalmip_sos_tools_missing";
    row.solve_time_sec = toc(t_case);
    cert.status = row.certificate_status;
    return;
end
if exist('linprog', 'file') ~= 2
    row.certificate_status = "solver_failed";
    row.solver_info = "linprog_missing";
    row.solve_time_sec = toc(t_case);
    cert.status = row.certificate_status;
    return;
end
if mod(dV, 2) ~= 0 || mod(dmult, 2) ~= 0
    row.certificate_status = "solver_failed";
    row.solver_info = "dV_and_dmult_must_be_even";
    row.solve_time_sec = toc(t_case);
    cert.status = row.certificate_status;
    return;
end
degree_meta = build_degree_meta(dV, dmult);
row = update_row_degree_meta(row, degree_meta);
if degree_meta.residual_has_odd_top_degree
    warning('Rev3:OddResidualPolynomialDegree', ...
        ['dV=%d,d_mult=%d gives residual_polynomial_degree=%d (odd). ', ...
        'The residual SOS Gram basis is explicitly built at d_sigma=%d, ', ...
        'but nonzero odd-degree residual terms must still be eliminated by the SDP.'], ...
        dV, dmult, degree_meta.residual_polynomial_degree, degree_meta.d_sigma);
end

reg_poly = scale_ctx.reg_poly;
[Phi, q] = build_sample_maps(reg_law, Ts, cfg.state_scale);
origin_eq = check_origin_equilibrium(q, scale_ctx.origin_regions, cfg.origin_equilibrium_tol);
forward_inv = verify_forward_invariance_union_rev3(reg_poly, Phi, q, scale_ctx, cfg);
row = update_row_forward_invariance(row, forward_inv);
row.partition_full_space_covered = scale_ctx.partition_full_space_covered;
row.full_space_partition_verified = scale_ctx.partition_full_space_covered;
row.origin_region_count = numel(scale_ctx.origin_regions);
row.origin_regions = join(string(scale_ctx.origin_regions), ";");
row = update_row_origin_equilibrium(row, origin_eq);

if ~(forward_inv.verified || forward_inv.tolerance_verified)
    if forward_inv.failed
        row.certificate_status = "forward_invariance_failed";
        row.solver_info = "skipped_due_to_forward_invariance_failure";
    else
        row.certificate_status = "forward_invariance_not_verified";
        row.solver_info = "skipped_due_to_forward_invariance_not_verified";
    end
    row.solve_time_sec = toc(t_case);
    cert.status = row.certificate_status;
    cert.row = row;
    cert.forward_invariance = forward_inv;
    cert.origin_regions = scale_ctx.origin_regions;
    cert.origin_equilibrium = origin_eq;
    return;
end
[transitions, trans_stats] = enumerate_transitions(reg_poly, Phi, q, cfg);
[partition_miss_count, partition_multi_count] = partition_grid_diagnostic(reg_poly, cfg.grid_x1, cfg.grid_x2);

size_meta = estimate_sos_size(reg_poly, transitions, degree_meta, dmult);
row = update_row_pre_solve(row, scale_ctx, transitions, trans_stats, ...
    partition_miss_count, partition_multi_count, degree_meta, size_meta);

if size_meta.sos_constraint_count > cfg.max_sos_constraints || ...
        size_meta.total_gram_dim_estimate > cfg.max_total_gram_dim
    row.certificate_status = "not_attempted_due_to_size";
    row.solver_info = "size_guard_triggered";
    row.solve_time_sec = toc(t_case);
    cert.status = row.certificate_status;
    cert.transitions = transitions;
    cert.size_meta = size_meta;
    return;
end

yalmip('clear');
tic;
[F, params, poly_meta] = build_sos_problem_rev3(reg_poly, Phi, q, transitions, ...
    scale_ctx.origin_regions, dV, dmult, degree_meta, cfg);
setup_time = toc;

ops = make_solver_options(cfg);
try
    [sol, ~, sos_grams, sos_residuals] = solvesos( ...
        F, [], ops, params, poly_meta.sos_candidate_monomials);
    sos_solver_return_block_count = numel(sos_grams);
catch ME
    sol = struct('problem', 999, 'info', string(getReport(ME, 'extended', 'hyperlinks', 'off')));
    sos_solver_return_block_count = NaN;
    sos_residuals = [];
end

row.setup_time_sec = setup_time;
row.solve_time_sec = toc(t_case);
row.solver_problem_code = sol.problem;
row.solver_info = string(sol.info);
row.yalmip_feasible = (sol.problem == 0);
sanity = empty_sanity();
dense = empty_dense_postcheck();
sos_oracle = empty_sos_oracle_postcheck();
origin_struct = empty_origin_polynomial_structure();

if sol.problem == 0
    [min_primal, postcheck_pass] = postcheck_constraints(F, cfg.postcheck_tol);
    t_sos_oracle = tic;
    sos_oracle = independent_sos_postcheck_rev3(poly_meta, cfg, sos_solver_return_block_count);
    sos_oracle.elapsed_time_sec = toc(t_sos_oracle);
    origin_struct = check_origin_polynomial_structure(poly_meta, scale_ctx.origin_regions, cfg.origin_structure_tol);
    sanity = sample_sanity_check(poly_meta, reg_poly, Phi, q, transitions, cfg);
    if cfg.dense_postcheck_enabled && Ts <= cfg.dense_postcheck_max_Ts + 10*eps
        t_dense = tic;
        dense = dense_grid_postcheck(poly_meta, reg_poly, Phi, q, transitions, cfg);
        dense.elapsed_time_sec = toc(t_dense);
    end
    if dense.attempted
        row.postcheck_method = "yalmip_check_plus_representative_point_sanity_plus_dense_grid";
    else
        row.postcheck_method = "yalmip_check_plus_representative_point_sanity";
    end
    row.postcheck_min_primal_residual = min_primal;
    row.postcheck_pass = postcheck_pass;
    row.sos_postcheck_method = sos_oracle.method;
    row.sos_postcheck_time_sec = sos_oracle.elapsed_time_sec;
    row.sos_postcheck_gram_source = sos_oracle.gram_source;
    row.sos_postcheck_gram_handle_missing_count = sos_oracle.gram_handle_missing_count;
    row.sos_postcheck_gram_dim_mismatch_count = sos_oracle.gram_dim_mismatch_count;
    row.sos_solver_return_block_count = sos_oracle.solver_return_block_count;
    row.gram_min_eig = sos_oracle.min_gram_eig;
    row.gram_worst_block = sos_oracle.worst_gram_block;
    row.coeff_match_max_abs = sos_oracle.max_abs_coeff_residual;
    row.coeff_match_worst_block = sos_oracle.worst_coeff_block;
    row.sos_oracle_strict_pass = sos_oracle.strict_pass;
    row.sos_oracle_tol_pass = sos_oracle.tol_pass;
    row = update_row_origin_polynomial_structure(row, origin_struct);
    row.sanity_min_region_pos = sanity.min_region_pos;
    row.sanity_min_transition_dec = sanity.min_transition_dec;
    row.sanity_worst_region = sanity.worst_region;
    row.sanity_worst_region_x1 = sanity.worst_region_x(1);
    row.sanity_worst_region_x2 = sanity.worst_region_x(2);
    row.sanity_worst_transition_index = sanity.worst_transition_index;
    row.sanity_worst_transition_i = sanity.worst_transition_i;
    row.sanity_worst_transition_j = sanity.worst_transition_j;
    row.sanity_worst_transition_x1 = sanity.worst_transition_x(1);
    row.sanity_worst_transition_x2 = sanity.worst_transition_x(2);
    row.sanity_worst_transition_xp1 = sanity.worst_transition_xp(1);
    row.sanity_worst_transition_xp2 = sanity.worst_transition_xp(2);
    row.sanity_pass = sanity.pass;
    row.sanity_strict_pass = sanity.strict_pass;
    row.sanity_tol_pass = sanity.tol_pass;
    row.dense_postcheck_attempted = dense.attempted;
    row.dense_postcheck_time_sec = dense.elapsed_time_sec;
    row.dense_grid_point_count = dense.grid_point_count;
    row.dense_min_region_pos = dense.min_region_pos;
    row.dense_min_transition_dec = dense.min_transition_dec;
    row.dense_worst_region = dense.worst_region;
    row.dense_worst_region_x1 = dense.worst_region_x(1);
    row.dense_worst_region_x2 = dense.worst_region_x(2);
    row.dense_worst_transition_index = dense.worst_transition_index;
    row.dense_worst_transition_i = dense.worst_transition_i;
    row.dense_worst_transition_j = dense.worst_transition_j;
    row.dense_worst_transition_x1 = dense.worst_transition_x(1);
    row.dense_worst_transition_x2 = dense.worst_transition_x(2);
    row.dense_worst_transition_xp1 = dense.worst_transition_xp(1);
    row.dense_worst_transition_xp2 = dense.worst_transition_xp(2);
    row.dense_strict_pass = dense.strict_pass;
    row.dense_tol_pass = dense.tol_pass;

    strict_postcheck_pass = sanity.strict_pass && (~dense.attempted || dense.strict_pass);
    tol_postcheck_pass = sanity.tol_pass && (~dense.attempted || dense.tol_pass);
    if forward_inv.tolerance_verified && ~forward_inv.verified
        row.certificate_status = "inconclusive";
    elseif postcheck_pass && strict_postcheck_pass && sos_oracle.strict_pass && ...
            forward_inv.verified && origin_eq.pass && origin_struct.pass
        row.certificate_status = "certified_asymptotic";
    elseif min_primal >= -10*cfg.postcheck_tol && tol_postcheck_pass && ...
            sos_oracle.tol_pass && forward_inv.verified && origin_eq.pass && origin_struct.pass
        row.certificate_status = "tol_certified_asymptotic";
    else
        row.certificate_status = "inconclusive";
    end
else
    if sol.problem == 1
        row.certificate_status = "inconclusive";
    else
        row.certificate_status = "solver_failed";
    end
end

cert.status = row.certificate_status;
cert.row = row;
cert.transitions = transitions;
cert.trans_stats = trans_stats;
cert.degree_meta = degree_meta;
cert.size_meta = size_meta;
cert.origin_regions = scale_ctx.origin_regions;
cert.origin_equilibrium = origin_eq;
cert.origin_polynomial_structure = origin_struct;
cert.forward_invariance = forward_inv;
cert.sanity = sanity;
cert.dense_postcheck = dense;
cert.sos_oracle = sos_oracle;
cert.sos_meta = summarize_sos_meta(poly_meta.sos_meta);
cert.sos_solver_residuals = sos_residuals;
cert.coefficients = extract_polynomial_coefficients(poly_meta);
end


function coeffs = extract_polynomial_coefficients(poly_meta)
M = numel(poly_meta.p_coeff);
coeffs = cell(1, M);
for i = 1:M
    coeffs{i} = value(poly_meta.p_coeff{i});
end
end


function sanity = sample_sanity_check(poly_meta, reg_poly, Phi, q, transitions, cfg)
x = poly_meta.x;
M = numel(reg_poly);
min_pos = inf;
min_dec = inf;
worst_region = NaN;
worst_region_x = [NaN; NaN];
worst_transition_index = NaN;
worst_transition_i = NaN;
worst_transition_j = NaN;
worst_transition_x = [NaN; NaN];
worst_transition_xp = [NaN; NaN];

for i = 1:M
    pts = representative_points(reg_poly(i).H, reg_poly(i).b, cfg);
    for pp = 1:size(pts, 2)
        xv = pts(:, pp);
        val = value(replace(poly_meta.p{i}, x, xv)) - cfg.epsilon * (xv' * xv);
        if val < min_pos
            min_pos = val;
            worst_region = i;
            worst_region_x = xv;
        end
    end
end

for kk = 1:numel(transitions)
    i = transitions(kk).i;
    j = transitions(kk).j;
    pts = representative_points(transitions(kk).G, transitions(kk).w, cfg);
    for pp = 1:size(pts, 2)
        xv = pts(:, pp);
        xp = Phi{i} * xv + q{i};
        val = value(replace(poly_meta.p{i}, x, xv)) - ...
            value(replace(poly_meta.p{j}, x, xp)) - cfg.alpha * (xv' * xv);
        if val < min_dec
            min_dec = val;
            worst_transition_index = kk;
            worst_transition_i = i;
            worst_transition_j = j;
            worst_transition_x = xv;
            worst_transition_xp = xp;
        end
    end
end

if isinf(min_dec)
    min_dec = NaN;
end
sanity = struct();
sanity.min_region_pos = min_pos;
sanity.min_transition_dec = min_dec;
sanity.worst_region = worst_region;
sanity.worst_region_x = worst_region_x;
sanity.worst_transition_index = worst_transition_index;
sanity.worst_transition_i = worst_transition_i;
sanity.worst_transition_j = worst_transition_j;
sanity.worst_transition_x = worst_transition_x;
sanity.worst_transition_xp = worst_transition_xp;
sanity.pass = min_pos >= -cfg.sanity_tol && (isnan(min_dec) || min_dec >= -cfg.sanity_tol);
sanity.tol_pass = min_pos >= -10*cfg.sanity_tol && (isnan(min_dec) || min_dec >= -10*cfg.sanity_tol);
sanity.strict_pass = min_pos >= 0 && (isnan(min_dec) || min_dec >= 0);
end


function sanity = empty_sanity()
sanity = struct();
sanity.min_region_pos = NaN;
sanity.min_transition_dec = NaN;
sanity.worst_region = NaN;
sanity.worst_region_x = [NaN; NaN];
sanity.worst_transition_index = NaN;
sanity.worst_transition_i = NaN;
sanity.worst_transition_j = NaN;
sanity.worst_transition_x = [NaN; NaN];
sanity.worst_transition_xp = [NaN; NaN];
sanity.pass = false;
sanity.tol_pass = false;
sanity.strict_pass = false;
end


function dense = dense_grid_postcheck(poly_meta, reg_poly, Phi, q, transitions, cfg)
x = poly_meta.x;
M = numel(reg_poly);
transitions_by_region = cell(1, M);
for kk = 1:numel(transitions)
    transitions_by_region{transitions(kk).i}(end+1) = kk; %#ok<AGROW>
end
min_pos = inf;
min_dec = inf;
worst_region = NaN;
worst_region_x = [NaN; NaN];
worst_transition_index = NaN;
worst_transition_i = NaN;
worst_transition_j = NaN;
worst_transition_x = [NaN; NaN];
worst_transition_xp = [NaN; NaN];
point_count = 0;

for a = 1:numel(cfg.dense_grid_scaled_x1)
    for b = 1:numel(cfg.dense_grid_scaled_x2)
        xv = [cfg.dense_grid_scaled_x1(a); cfg.dense_grid_scaled_x2(b)];
        point_count = point_count + 1;
        hit_regions = [];
        for i = 1:M
            if all(reg_poly(i).H * xv <= reg_poly(i).b + cfg.lp_tol)
                hit_regions(end+1) = i; %#ok<AGROW>
                val = value(replace(poly_meta.p{i}, x, xv)) - cfg.epsilon * (xv' * xv);
                if val < min_pos
                    min_pos = val;
                    worst_region = i;
                    worst_region_x = xv;
                end
            end
        end
        for hh = 1:numel(hit_regions)
            candidate_transitions = transitions_by_region{hit_regions(hh)};
            for cc = 1:numel(candidate_transitions)
                kk = candidate_transitions(cc);
                if ~all(transitions(kk).G * xv <= transitions(kk).w + cfg.lp_tol)
                    continue;
                end
                i = transitions(kk).i;
                j = transitions(kk).j;
                xp = Phi{i} * xv + q{i};
                val = value(replace(poly_meta.p{i}, x, xv)) - ...
                    value(replace(poly_meta.p{j}, x, xp)) - cfg.alpha * (xv' * xv);
                if val < min_dec
                    min_dec = val;
                    worst_transition_index = kk;
                    worst_transition_i = i;
                    worst_transition_j = j;
                    worst_transition_x = xv;
                    worst_transition_xp = xp;
                end
            end
        end
    end
end

if isinf(min_pos)
    min_pos = NaN;
end
if isinf(min_dec)
    min_dec = NaN;
end

dense = struct();
dense.attempted = true;
dense.elapsed_time_sec = NaN;
dense.grid_point_count = point_count;
dense.min_region_pos = min_pos;
dense.min_transition_dec = min_dec;
dense.worst_region = worst_region;
dense.worst_region_x = worst_region_x;
dense.worst_transition_index = worst_transition_index;
dense.worst_transition_i = worst_transition_i;
dense.worst_transition_j = worst_transition_j;
dense.worst_transition_x = worst_transition_x;
dense.worst_transition_xp = worst_transition_xp;
dense.strict_pass = min_pos >= 0 && (isnan(min_dec) || min_dec >= 0);
dense.tol_pass = min_pos >= -10*cfg.sanity_tol && (isnan(min_dec) || min_dec >= -10*cfg.sanity_tol);
end


function dense = empty_dense_postcheck()
dense = struct();
dense.attempted = false;
dense.elapsed_time_sec = NaN;
dense.grid_point_count = 0;
dense.min_region_pos = NaN;
dense.min_transition_dec = NaN;
dense.worst_region = NaN;
dense.worst_region_x = [NaN; NaN];
dense.worst_transition_index = NaN;
dense.worst_transition_i = NaN;
dense.worst_transition_j = NaN;
dense.worst_transition_x = [NaN; NaN];
dense.worst_transition_xp = [NaN; NaN];
dense.strict_pass = false;
dense.tol_pass = false;
end


function pts = representative_points(G, w, cfg)
pts = zeros(2, 0);
[is_feas, x_feas] = poly_feasible(G, w, cfg.lp_tol);
if is_feas
    pts(:, end+1) = x_feas; %#ok<AGROW>
end
verts = enumerate_vertices(G, w, cfg.lp_tol);
if ~isempty(verts)
    pts = [pts, verts]; %#ok<AGROW>
end
if isempty(pts)
    pts = zeros(2, 1);
end
pts = unique(round(pts', 10), 'rows')';
end


function [min_primal, pass] = postcheck_constraints(F, tol)
try
    primal = check(F);
    if isempty(primal)
        min_primal = NaN;
        pass = false;
    else
        min_primal = min(primal);
        pass = all(primal >= -tol);
    end
catch
    min_primal = NaN;
    pass = false;
end
end


function sos_oracle = empty_sos_oracle_postcheck()
sos_oracle = struct();
sos_oracle.method = "not_run";
sos_oracle.elapsed_time_sec = NaN;
sos_oracle.min_gram_eig = NaN;
sos_oracle.max_abs_coeff_residual = NaN;
sos_oracle.worst_gram_block = "";
sos_oracle.worst_coeff_block = "";
sos_oracle.sos_block_count = 0;
sos_oracle.gram_source = "";
sos_oracle.gram_handle_missing_count = 0;
sos_oracle.gram_dim_mismatch_count = 0;
sos_oracle.solver_return_block_count = NaN;
sos_oracle.strict_pass = false;
sos_oracle.tol_pass = false;
end


function summary = summarize_sos_meta(sos_meta)
summary = struct('name', {}, 'kind', {}, 'degree', {}, 'gram_dim', {}, ...
    'has_gram_handle', {}, 'gram_source', {}, 'gram_monomial_count', {}, ...
    'candidate_monomial_count', {}, 'gram_basis_mismatch', {});
for kk = 1:numel(sos_meta)
    summary(kk).name = sos_meta(kk).name; %#ok<AGROW>
    summary(kk).kind = sos_meta(kk).kind;
    summary(kk).degree = sos_meta(kk).degree;
    summary(kk).gram_dim = sos_meta(kk).gram_dim;
    summary(kk).has_gram_handle = ~isempty(sos_meta(kk).sos_constraint);
    summary(kk).gram_source = sos_meta(kk).gram_source;
    summary(kk).gram_monomial_count = sos_meta(kk).candidate_monomial_count;
    summary(kk).candidate_monomial_count = sos_meta(kk).candidate_monomial_count;
    summary(kk).gram_basis_mismatch = sos_meta(kk).gram_basis_mismatch;
end
end


function degree_meta = build_degree_meta(dV, dmult)
degree_meta = struct();
degree_meta.d_V = dV;
degree_meta.d_mult = dmult;
degree_meta.residual_polynomial_degree = max(dV, dmult + 1);
degree_meta.d_sigma = smallest_even_at_least(max(dV, dmult + 1));
degree_meta.residual_sos_degree = degree_meta.d_sigma;
degree_meta.residual_has_odd_top_degree = mod(degree_meta.residual_polynomial_degree, 2) ~= 0;
degree_meta.lyapunov_monomial_count = monomial_count_2d(dV);
degree_meta.residual_gram_dim = monomial_count_2d(degree_meta.d_sigma / 2);
if dmult == 0
    degree_meta.multiplier_gram_dim = 1;
else
    degree_meta.multiplier_gram_dim = monomial_count_2d(dmult / 2);
end
end


function size_meta = estimate_sos_size(reg_poly, transitions, degree_meta, dmult)
region_generator_count = 0;
for i = 1:numel(reg_poly)
    region_generator_count = region_generator_count + size(reg_poly(i).H, 1);
end
transition_generator_count = 0;
for kk = 1:numel(transitions)
    transition_generator_count = transition_generator_count + size(transitions(kk).G, 1);
end
residual_count = numel(reg_poly) + numel(transitions);
multiplier_sos_count = region_generator_count + transition_generator_count;
size_meta = struct();
size_meta.region_generator_count = region_generator_count;
size_meta.transition_generator_count = transition_generator_count;
size_meta.residual_sos_count = residual_count;
size_meta.multiplier_sos_count = multiplier_sos_count;
size_meta.sos_constraint_count = residual_count + multiplier_sos_count;
size_meta.total_gram_dim_estimate = residual_count * degree_meta.residual_gram_dim + ...
    multiplier_sos_count * degree_meta.multiplier_gram_dim;
end


function d = smallest_even_at_least(v)
d = ceil(v);
if mod(d, 2) ~= 0
    d = d + 1;
end
end


function n = monomial_count_2d(d)
n = nchoosek(2 + d, d);
end


function row = base_row(Ts, dV, dmult, cfg)
row = struct();
row.Ts = Ts;
row.d_V = dV;
row.d_mult = dmult;
row.d_sigma = smallest_even_at_least(max(dV, dmult + 1));
row.residual_polynomial_degree = max(dV, dmult + 1);
row.residual_sos_degree = row.d_sigma;
row.residual_has_odd_top_degree = mod(row.residual_polynomial_degree, 2) ~= 0;
row.epsilon = cfg.epsilon;
row.alpha = cfg.alpha;
row.solver = string(cfg.solver);
row.run_label = string(cfg.run_label);
row.use_parallel = cfg.use_parallel;
row.parallel_workers = cfg.parallel_workers;
row.yalmip_sos_model = cfg.yalmip_sos_model;
row.yalmip_sos_scale = cfg.yalmip_sos_scale;
row.yalmip_sos_newton = cfg.yalmip_sos_newton;
row.mosek_primal_feas_tol = cfg.mosek_primal_feas_tol;
row.mosek_dual_feas_tol = cfg.mosek_dual_feas_tol;
row.mosek_rel_gap_tol = cfg.mosek_rel_gap_tol;
row.mosek_max_iter = cfg.mosek_max_iter;
row.mosek_num_threads = cfg.mosek_num_threads;
row.state_scale_s1 = cfg.state_scale(1);
row.state_scale_s2 = cfg.state_scale(2);
row.normalization_mode = string(cfg.normalization_mode);
row.coefficient_bound = cfg.coefficient_bound;
row.certificate_type = "asymptotic";
row.certificate_status = "inconclusive";
row.solver_problem_code = NaN;
row.solver_info = "";
row.yalmip_feasible = false;
row.postcheck_method = "not_run";
row.postcheck_min_primal_residual = NaN;
row.postcheck_pass = false;
row.sos_postcheck_method = "not_run";
row.sos_postcheck_time_sec = NaN;
row.sos_postcheck_gram_source = "";
row.sos_postcheck_gram_handle_missing_count = NaN;
row.sos_postcheck_gram_dim_mismatch_count = NaN;
row.sos_solver_return_block_count = NaN;
row.gram_min_eig = NaN;
row.gram_worst_block = "";
row.coeff_match_max_abs = NaN;
row.coeff_match_worst_block = "";
row.sos_oracle_strict_pass = false;
row.sos_oracle_tol_pass = false;
row.sanity_min_region_pos = NaN;
row.sanity_min_transition_dec = NaN;
row.sanity_worst_region = NaN;
row.sanity_worst_region_x1 = NaN;
row.sanity_worst_region_x2 = NaN;
row.sanity_worst_transition_index = NaN;
row.sanity_worst_transition_i = NaN;
row.sanity_worst_transition_j = NaN;
row.sanity_worst_transition_x1 = NaN;
row.sanity_worst_transition_x2 = NaN;
row.sanity_worst_transition_xp1 = NaN;
row.sanity_worst_transition_xp2 = NaN;
row.sanity_pass = false;
row.sanity_strict_pass = false;
row.sanity_tol_pass = false;
row.dense_grid_coordinate_system = string(cfg.dense_grid_coordinate_system);
row.dense_grid_scaled_x1_min = min(cfg.dense_grid_scaled_x1);
row.dense_grid_scaled_x1_max = max(cfg.dense_grid_scaled_x1);
row.dense_grid_scaled_x1_count = numel(cfg.dense_grid_scaled_x1);
row.dense_grid_scaled_x2_min = min(cfg.dense_grid_scaled_x2);
row.dense_grid_scaled_x2_max = max(cfg.dense_grid_scaled_x2);
row.dense_grid_scaled_x2_count = numel(cfg.dense_grid_scaled_x2);
row.dense_postcheck_attempted = false;
row.dense_postcheck_time_sec = NaN;
row.dense_grid_point_count = 0;
row.dense_min_region_pos = NaN;
row.dense_min_transition_dec = NaN;
row.dense_worst_region = NaN;
row.dense_worst_region_x1 = NaN;
row.dense_worst_region_x2 = NaN;
row.dense_worst_transition_index = NaN;
row.dense_worst_transition_i = NaN;
row.dense_worst_transition_j = NaN;
row.dense_worst_transition_x1 = NaN;
row.dense_worst_transition_x2 = NaN;
row.dense_worst_transition_xp1 = NaN;
row.dense_worst_transition_xp2 = NaN;
row.dense_strict_pass = false;
row.dense_tol_pass = false;
row.forward_invariance_status = "";
row.forward_invariance_verified = false;
row.forward_invariance_tolerance_verified = false;
row.forward_invariance_checked = false;
row.forward_invariance_method = "";
row.forward_invariance_lp_count = NaN;
row.forward_invariance_max_escape_margin = NaN;
row.forward_invariance_max_escape_margin_region = NaN;
row.forward_invariance_max_escape_margin_choice = "";
row.forward_invariance_escape_margin_tol = cfg.forward_invariance_tol;
row.forward_invariance_counterexample_region = NaN;
row.forward_invariance_counterexample_x1 = NaN;
row.forward_invariance_counterexample_x2 = NaN;
row.forward_invariance_counterexample_xp1 = NaN;
row.forward_invariance_counterexample_xp2 = NaN;
row.forward_invariance_counterexample_choice = "";
row.partition_full_space_covered = false;
row.full_space_partition_verified = false;
row.origin_region_count = NaN;
row.origin_regions = "";
row.origin_equilibrium_status = "not_checked";
row.max_origin_q_norm = NaN;
row.origin_structure_status = "not_checked";
row.max_origin_p_abs = NaN;
row.max_origin_grad_norm = NaN;
row.reachable_transition_count = NaN;
row.ambiguous_transition_count = NaN;
row.degenerate_transition_count = NaN;
row.near_degenerate_transition_count = NaN;
row.min_transition_cheby_radius = NaN;
row.partition_grid_miss_count = NaN;
row.partition_grid_multi_count = NaN;
row.lyapunov_monomial_count = NaN;
row.residual_gram_dim = NaN;
row.multiplier_gram_dim = NaN;
row.sos_constraint_count = NaN;
row.total_gram_dim_estimate = NaN;
row.region_generator_count = NaN;
row.transition_generator_count = NaN;
row.setup_time_sec = NaN;
row.solve_time_sec = NaN;
end


function row = update_row_degree_meta(row, degree_meta)
row.d_sigma = degree_meta.d_sigma;
row.residual_polynomial_degree = degree_meta.residual_polynomial_degree;
row.residual_sos_degree = degree_meta.residual_sos_degree;
row.residual_has_odd_top_degree = degree_meta.residual_has_odd_top_degree;
row.lyapunov_monomial_count = degree_meta.lyapunov_monomial_count;
row.residual_gram_dim = degree_meta.residual_gram_dim;
row.multiplier_gram_dim = degree_meta.multiplier_gram_dim;
end


function row = update_row_pre_solve(row, scale_ctx, transitions, trans_stats, ...
    partition_miss_count, partition_multi_count, degree_meta, size_meta)
row.partition_full_space_covered = scale_ctx.partition_full_space_covered;
row.full_space_partition_verified = scale_ctx.partition_full_space_covered;
row.origin_region_count = numel(scale_ctx.origin_regions);
row.origin_regions = join(string(scale_ctx.origin_regions), ";");
row.reachable_transition_count = numel(transitions);
row.ambiguous_transition_count = trans_stats.ambiguous_count;
row.degenerate_transition_count = trans_stats.degenerate_count;
row.near_degenerate_transition_count = trans_stats.near_degenerate_count;
row.min_transition_cheby_radius = trans_stats.min_radius;
row.partition_grid_miss_count = partition_miss_count;
row.partition_grid_multi_count = partition_multi_count;
row = update_row_degree_meta(row, degree_meta);
row.sos_constraint_count = size_meta.sos_constraint_count;
row.total_gram_dim_estimate = size_meta.total_gram_dim_estimate;
row.region_generator_count = size_meta.region_generator_count;
row.transition_generator_count = size_meta.transition_generator_count;
end


function row = update_row_forward_invariance(row, forward_inv)
row.forward_invariance_status = forward_inv.status;
row.forward_invariance_verified = forward_inv.verified;
row.forward_invariance_tolerance_verified = forward_inv.tolerance_verified;
row.forward_invariance_checked = forward_inv.checked;
row.forward_invariance_method = forward_inv.method;
row.forward_invariance_lp_count = forward_inv.lp_count;
row.forward_invariance_max_escape_margin = forward_inv.max_escape_margin;
row.forward_invariance_max_escape_margin_region = forward_inv.max_escape_margin_region;
row.forward_invariance_max_escape_margin_choice = forward_inv.max_escape_margin_choice;
row.forward_invariance_counterexample_region = forward_inv.counterexample_region;
row.forward_invariance_counterexample_x1 = forward_inv.counterexample_x(1);
row.forward_invariance_counterexample_x2 = forward_inv.counterexample_x(2);
row.forward_invariance_counterexample_xp1 = forward_inv.counterexample_xp(1);
row.forward_invariance_counterexample_xp2 = forward_inv.counterexample_xp(2);
row.forward_invariance_counterexample_choice = forward_inv.counterexample_choice;
end


function row = update_row_origin_equilibrium(row, origin_eq)
row.origin_equilibrium_status = origin_eq.status;
row.max_origin_q_norm = origin_eq.max_q_norm;
end


function row = update_row_origin_polynomial_structure(row, origin_struct)
row.origin_structure_status = origin_struct.status;
row.max_origin_p_abs = origin_struct.max_p_abs;
row.max_origin_grad_norm = origin_struct.max_grad_norm;
end


function cert = empty_cert(Ts, dV, dmult, cfg)
cert = struct();
cert.Ts = Ts;
cert.d_V = dV;
cert.d_mult = dmult;
cert.epsilon = cfg.epsilon;
cert.alpha = cfg.alpha;
cert.status = "not_started";
end


function rows = append_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end


function scale_ctx = build_scale_context(reg, cfg)
T = diag(cfg.state_scale);
scale_ctx = struct();
scale_ctx.state_scale_matrix = T;
[scale_ctx.reg_poly, scale_ctx.region_zero_row_removed_count, ...
    scale_ctx.removed_region_constraint_count] = scale_region_polytopes(reg, T, cfg);
[~, scale_ctx.origin_regions] = find_origin_regions(scale_ctx.reg_poly, cfg.lp_tol);
scale_ctx.partition_proof = verify_full_state_space_partition(scale_ctx.reg_poly, cfg.partition_tol);
scale_ctx.partition_full_space_covered = scale_ctx.partition_proof.full_space_covered;
end


function [reg_poly, zero_removed_total, removed_region_total] = scale_region_polytopes(reg, T, cfg)
reg_poly = reg;
zero_removed_total = 0;
removed_region_total = 0;
for i = 1:numel(reg)
    H = reg(i).H * T;
    b = reg(i).b;
    [H, b, zero_removed, removed_redundant, is_empty] = preprocess_polyhedron(H, b, cfg);
    if is_empty
        error('Region %d became empty during preprocessing.', i);
    end
    reg_poly(i).H = H;
    reg_poly(i).b = b;
    zero_removed_total = zero_removed_total + zero_removed;
    removed_region_total = removed_region_total + removed_redundant;
end
end


function [Phi, q] = build_sample_maps(reg, Ts, state_scale)
M = numel(reg);
Phi = cell(1, M);
q = cell(1, M);
A0 = [1 Ts; 0 1];
B0 = [Ts^2/2; Ts];
T = diag(state_scale);

for i = 1:M
    Phi_original = A0 + B0 * reg(i).K;
    q_original = B0 * reg(i).k;
    Phi{i} = T \ (Phi_original * T);
    q{i} = T \ q_original;
end
end


function origin_eq = check_origin_equilibrium(q, origin_regions, tol)
origin_eq = struct();
origin_eq.origin_regions = origin_regions;
origin_eq.q_norms = [];
origin_eq.max_q_norm = NaN;
origin_eq.pass = false;
origin_eq.status = "not_checked";

if isempty(origin_regions)
    origin_eq.status = "no_origin_region";
    return;
end

q_norms = zeros(1, numel(origin_regions));
for kk = 1:numel(origin_regions)
    i = origin_regions(kk);
    q_norms(kk) = norm(q{i}, 2);
end

origin_eq.q_norms = q_norms;
origin_eq.max_q_norm = max(q_norms);
origin_eq.pass = origin_eq.max_q_norm <= tol;
if origin_eq.pass
    origin_eq.status = "verified_origin_equilibrium";
else
    origin_eq.status = "failed_origin_equilibrium";
end
end


function origin_struct = empty_origin_polynomial_structure()
origin_struct = struct();
origin_struct.origin_regions = [];
origin_struct.p_abs = [];
origin_struct.grad_norms = [];
origin_struct.max_p_abs = NaN;
origin_struct.max_grad_norm = NaN;
origin_struct.pass = false;
origin_struct.status = "not_checked";
end


function origin_struct = check_origin_polynomial_structure(poly_meta, origin_regions, tol)
origin_struct = empty_origin_polynomial_structure();
origin_struct.origin_regions = origin_regions;

if isempty(origin_regions)
    origin_struct.status = "no_origin_region";
    return;
end

x = poly_meta.x;
x0 = zeros(numel(x), 1);
p_abs = zeros(1, numel(origin_regions));
grad_norms = zeros(1, numel(origin_regions));
for kk = 1:numel(origin_regions)
    i = origin_regions(kk);
    p0 = value(replace(poly_meta.p{i}, x, x0));
    gradp = jacobian(poly_meta.p{i}, x);
    grad0 = value(replace(gradp, x, x0));
    p_abs(kk) = abs(p0);
    grad_norms(kk) = norm(grad0, 2);
end

origin_struct.p_abs = p_abs;
origin_struct.grad_norms = grad_norms;
origin_struct.max_p_abs = max(p_abs);
origin_struct.max_grad_norm = max(grad_norms);
origin_struct.pass = all(isfinite(p_abs)) && all(isfinite(grad_norms)) && ...
    origin_struct.max_p_abs <= tol && origin_struct.max_grad_norm <= tol;
if origin_struct.pass
    origin_struct.status = "verified_origin_polynomial_structure";
else
    origin_struct.status = "failed_origin_polynomial_structure";
end
end


function [Gout, wout, zero_removed, removed_redundant, is_empty] = preprocess_polyhedron(G, w, cfg)
zero_rows = vecnorm(G, 2, 2) <= cfg.zero_row_tol;
if any(zero_rows & (w < -cfg.lp_tol))
    Gout = G;
    wout = w;
    zero_removed = 0;
    removed_redundant = 0;
    is_empty = true;
    return;
end
keep = ~zero_rows;
G = G(keep, :);
w = w(keep);
zero_removed = nnz(zero_rows);
[is_feas] = poly_feasible(G, w, cfg.lp_tol);
if ~is_feas
    Gout = G;
    wout = w;
    removed_redundant = 0;
    is_empty = true;
    return;
end
[Gout, wout, removed_redundant] = remove_redundant_constraints(G, w, cfg.redundancy_tol);
is_empty = false;
end


function [transitions, stats] = enumerate_transitions(reg, Phi, q, cfg)
transitions = struct('i', {}, 'j', {}, 'G', {}, 'w', {}, 'cheby_radius', {}, 'x_feas', {}, 'is_ambiguous', {});
stats = struct('min_radius', inf, 'degenerate_count', 0, ...
    'near_degenerate_count', 0, 'ambiguous_count', 0);

for i = 1:numel(reg)
    for j = 1:numel(reg)
        Gij = [reg(i).H; reg(j).H * Phi{i}];
        wij = [reg(i).b; reg(j).b - reg(j).H * q{i}];
        [Gij, wij, ~, ~, is_empty] = preprocess_polyhedron(Gij, wij, cfg);
        if is_empty
            continue;
        end
        [is_feas, x_feas, radius, is_ambiguous] = poly_feasibility_and_cheby_radius(Gij, wij, cfg.lp_tol);
        if is_feas || is_ambiguous
            transitions(end+1) = struct('i', i, 'j', j, 'G', Gij, 'w', wij, ...
                'cheby_radius', radius, 'x_feas', x_feas, 'is_ambiguous', is_ambiguous); %#ok<AGROW>
            if isfinite(radius)
                stats.min_radius = min(stats.min_radius, radius);
                if radius <= cfg.lp_tol
                    stats.degenerate_count = stats.degenerate_count + 1;
                elseif radius <= cfg.near_degenerate_radius_tol
                    stats.near_degenerate_count = stats.near_degenerate_count + 1;
                end
            end
            if is_ambiguous
                stats.ambiguous_count = stats.ambiguous_count + 1;
            end
        end
    end
end
if isinf(stats.min_radius)
    stats.min_radius = NaN;
end
end


function [is_feas, x_feas, radius, is_ambiguous] = poly_feasibility_and_cheby_radius(G, w, tol)
n = size(G, 2);
x_feas = NaN(n, 1);
radius = NaN;
is_ambiguous = false;

if isempty(G)
    is_feas = true;
    x_feas = zeros(n, 1);
    radius = inf;
    return;
end

norms = vecnorm(G, 2, 2);
A = [G, norms];
b = w;
f = [zeros(n, 1); -1];
lb = [-inf(n, 1); 0];
opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
[sol, fval, exitflag] = linprog(f, A, b, [], [], lb, [], opts);
if exitflag == 1
    is_feas = true;
    x_feas = sol(1:n);
    radius = -fval;
    return;
end

[x0, ~, exitflag0] = linprog(zeros(n,1), G, w, [], [], [], [], opts);
if exitflag0 == 1
    is_feas = true;
    x_feas = x0;
    radius = 0;
    return;
end

[x2, ~, exitflag2] = linprog(zeros(n,1), G, w + tol, [], [], [], [], opts);
if exitflag2 == 1
    is_feas = true;
    x_feas = x2;
    radius = 0;
    is_ambiguous = true;
else
    is_feas = false;
end
end


function [is_feas, x_feas] = poly_feasible(G, w, tol)
if isempty(G)
    is_feas = true;
    x_feas = zeros(2, 1);
    return;
end
opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
[x_feas, ~, exitflag] = linprog(zeros(size(G,2),1), G, w + tol, [], [], [], [], opts);
is_feas = (exitflag == 1);
if ~is_feas
    x_feas = NaN(size(G,2), 1);
end
end


function [Gkeep, wkeep, removed_count] = remove_redundant_constraints(G, w, tol)
m = size(G, 1);
keep = true(m, 1);
opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
for r = 1:m
    other = keep;
    other(r) = false;
    if ~any(other)
        continue;
    end
    f = -G(r, :)';
    [~, fval, exitflag] = linprog(f, G(other, :), w(other), [], [], [], [], opts);
    if exitflag == 1
        max_val = -fval;
        if max_val <= w(r) + tol
            keep(r) = false;
        end
    end
end
Gkeep = G(keep, :);
wkeep = w(keep);
removed_count = nnz(~keep);
end


function [origin_mask, origin_regions] = find_origin_regions(reg_poly, tol)
M = numel(reg_poly);
origin_mask = false(1, M);
for i = 1:M
    origin_mask(i) = all(reg_poly(i).H * [0; 0] <= reg_poly(i).b + tol);
end
origin_regions = find(origin_mask);
end


function proof = verify_full_state_space_partition(reg, tol)
M = numel(reg);
choices = ones(1, M);
best_x = [NaN; NaN];
best_choice = choices;
best_margin = -inf;
checked = 0;
invalid_candidate_count = 0;
ambiguous_lp_count = 0;
opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
while true
    G = zeros(M, 2);
    w = zeros(M, 1);
    for i = 1:M
        row = choices(i);
        G(i, :) = -reg(i).H(row, :);
        w(i) = -reg(i).b(row) - tol;
    end
    [x, ~, exitflag] = linprog([0; 0], G, w, [], [], [], [], opts);
    checked = checked + 1;
    if exitflag == 1 || exitflag == -3
        [valid_counterexample, complement_margin] = ...
            is_valid_full_space_counterexample(reg, x, tol);
        if isfinite(complement_margin) && complement_margin > best_margin
            best_margin = complement_margin;
            best_x = x;
            best_choice = choices;
        end
        if ~valid_counterexample
            invalid_candidate_count = invalid_candidate_count + 1;
            [choices, done] = next_full_space_choice_vector(choices, reg);
            if done
                break;
            end
            continue;
        end
        proof = struct('full_space_covered', false, 'checked_cases', checked, ...
            'counterexample', x, 'max_complement_margin', complement_margin, ...
            'best_x', best_x, 'best_choice', best_choice, ...
            'invalid_candidate_count', invalid_candidate_count, ...
            'ambiguous_lp_count', ambiguous_lp_count);
        return;
    elseif exitflag ~= -2
        ambiguous_lp_count = ambiguous_lp_count + 1;
    end
    [choices, done] = next_full_space_choice_vector(choices, reg);
    if done
        break;
    end
end
proof = struct('full_space_covered', ambiguous_lp_count == 0, ...
    'checked_cases', checked, ...
    'counterexample', best_x, 'max_complement_margin', best_margin, ...
    'best_x', best_x, 'best_choice', best_choice, ...
    'invalid_candidate_count', invalid_candidate_count, ...
    'ambiguous_lp_count', ambiguous_lp_count);
end


function [tf, complement_margin] = is_valid_full_space_counterexample(reg, x, tol)
tf = false;
complement_margin = NaN;
if isempty(x) || numel(x) ~= size(reg(1).H, 2) || any(~isfinite(x))
    return;
end
region_margins = zeros(1, numel(reg));
for i = 1:numel(reg)
    region_margins(i) = max(reg(i).H * x(:) - reg(i).b);
end
complement_margin = min(region_margins);
tf = complement_margin > tol;
end


function [choices, done] = next_full_space_choice_vector(choices, reg)
done = false;
k = 1;
while k <= numel(choices)
    choices(k) = choices(k) + 1;
    if choices(k) <= size(reg(k).H, 1)
        return;
    end
    choices(k) = 1;
    k = k + 1;
end
done = true;
end


function [miss_count, multi_count] = partition_grid_diagnostic(reg_poly, grid_x1, grid_x2)
miss_count = 0;
multi_count = 0;
for a = 1:numel(grid_x1)
    for b = 1:numel(grid_x2)
        x = [grid_x1(a); grid_x2(b)];
        hits = 0;
        for i = 1:numel(reg_poly)
            if all(reg_poly(i).H * x <= reg_poly(i).b + 1e-8)
                hits = hits + 1;
            end
        end
        if hits == 0
            miss_count = miss_count + 1;
        elseif hits > 1
            multi_count = multi_count + 1;
        end
    end
end
end


function vertices = enumerate_vertices(G, w, tol)
n = size(G, 2);
if n ~= 2 || size(G, 1) < 2
    vertices = zeros(2, 0);
    return;
end
vertices = zeros(2, 0);
for a = 1:size(G, 1)-1
    for b = a+1:size(G, 1)
        Aeq = [G(a, :); G(b, :)];
        if abs(det(Aeq)) <= tol
            continue;
        end
        x = Aeq \ [w(a); w(b)];
        if all(G*x <= w + 10*tol)
            vertices(:, end+1) = x; %#ok<AGROW>
        end
    end
end
if ~isempty(vertices)
    vertices = unique(round(vertices', 10), 'rows')';
end
end


function script_dir = get_script_dir()
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
end


function root = find_release_root(start_dir)
root = start_dir;
while true
    if exist(fullfile(root, 'data', 'double_integrator', 'answer_regions.m'), 'file') == 2
        return;
    end
    parent = fileparts(root);
    if strcmp(parent, root)
        error('Could not find release root from %s', start_dir);
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


function ensure_yalmip_export_precedence()
% MATLAB ships a codetools/export.m that can shadow YALMIP's extras/export.m.
% solvesos/compilesos needs the YALMIP export implementation.
sdpvar_path = which('sdpvar');
if isempty(sdpvar_path)
    return;
end
yalmip_root = fileparts(fileparts(sdpvar_path));
yalmip_extras = fullfile(yalmip_root, 'extras');
if exist(fullfile(yalmip_extras, 'export.m'), 'file') == 2
    addpath(yalmip_extras, '-begin');
end
end


function configure_mosek_path(cfg)
if ~isfield(cfg, 'mosek_path') || strlength(string(cfg.mosek_path)) == 0
    return;
end
mosek_path = char(cfg.mosek_path);
if isfolder(mosek_path) && ~any(strcmpi(strsplit(path, pathsep), mosek_path))
    addpath(mosek_path);
end
end


function ensure_matlab_linprog_precedence()
optim_path = fullfile(matlabroot, 'toolbox', 'optim', 'optim');
if isfolder(optim_path)
    addpath(optim_path, '-begin');
end
end


function ops = make_solver_options(cfg)
ops = sdpsettings('solver', char(cfg.solver), 'verbose', cfg.verbose);
ops.cachesolvers = 1;
ops.sos.model = cfg.yalmip_sos_model;
ops.sos.scale = cfg.yalmip_sos_scale;
ops.sos.newton = cfg.yalmip_sos_newton;

if strcmpi(string(cfg.solver), "mosek")
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = cfg.mosek_primal_feas_tol;
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = cfg.mosek_dual_feas_tol;
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = cfg.mosek_rel_gap_tol;
    ops.mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS = cfg.mosek_max_iter;
    ops.mosek.MSK_IPAR_NUM_THREADS = cfg.mosek_num_threads;
end
end
