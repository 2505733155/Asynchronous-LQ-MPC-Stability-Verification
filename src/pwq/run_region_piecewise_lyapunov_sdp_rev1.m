%% Rev1 fixed-Ts explicit-MPC PWA Lyapunov SDP scan
% Stable one-step regional/sectional quadratic Lyapunov SDP.
%
% Rev1 starts from an early one-step PWA/PWQ Lyapunov baseline and folds in
% only the useful numerical preprocessing found in the polytope-refinement
% experiments:
% state scaling, all-origin-region handling, conservative transition
% enumeration, redundancy removal, and post-solve LMI checks.

clear; clc;

cfg.sample_times = [0.80 0.82 0.83 0.84 0.85 0.86 0.88 0.90 0.95 1.00 1.25 1.50];
cfg.margin_values = [1e-9 1e-11];
cfg.multiplier_modes = ["scalar", "matrix"];
cfg.solver = "sdpt3";
cfg.verbose = 0;
cfg.lp_tol = 1e-8;
cfg.postcheck_tol = 1e-7;
cfg.grid_x1 = linspace(-30, 30, 121);
cfg.grid_x2 = linspace(-12, 12, 97);
cfg.sanity_tol = 1e-8;
cfg.dense_postcheck_enabled = true;
cfg.dense_postcheck_max_Ts = 1.0;

cfg.state_scale = [20, 5];
cfg.dense_grid_coordinate_system = "scaled_state_from_physical_box";
cfg.dense_grid_scaled_x1 = linspace(-100 / cfg.state_scale(1), 100 / cfg.state_scale(1), 15);
cfg.dense_grid_scaled_x2 = linspace(-100 / cfg.state_scale(2), 100 / cfg.state_scale(2), 57);
cfg.row_normalization_mode = "none";
cfg.redundancy_mode = "remove_redundant";
cfg.equality_multiplier_mode = "none";
cfg.boundary_continuity_mode = "off";
cfg.zero_row_tol = 1e-12;
cfg.redundancy_tol = 1e-9;
cfg.near_degenerate_radius_tol = 1e-3;
cfg.preprocess_safety_status = build_preprocess_safety_status(cfg);
cfg.mosek_path = getenv('AMPC_MOSEK_PATH');
cfg.mosek_primal_feas_tol = 1e-8;
cfg.mosek_dual_feas_tol = 1e-8;
cfg.mosek_rel_gap_tol = 1e-8;
cfg.mosek_max_iter = 200;
cfg.mosek_num_threads = 1;
cfg.result_filename = 'region_piecewise_lyapunov_sdp_rev1_results.csv';

if strcmpi(getenv('PWQ_REV1_FAIR_MOSEK'), '1')
    cfg.sample_times = 0.1:0.1:1.0;
    cfg.solver = "mosek";
    cfg.result_filename = 'region_piecewise_lyapunov_sdp_rev1_mosek_fair_results.csv';
end

shared_mosek_tol = getenv('PWQ_REV1_MOSEK_TOL');
if ~isempty(shared_mosek_tol)
    tol = str2double(shared_mosek_tol);
    if ~isfinite(tol) || tol <= 0
        error('Invalid PWQ_REV1_MOSEK_TOL: %s', shared_mosek_tol);
    end
    cfg.mosek_primal_feas_tol = tol;
    cfg.mosek_dual_feas_tol = tol;
    cfg.mosek_rel_gap_tol = tol;
end

result_filename_override = getenv('PWQ_REV1_RESULT_FILENAME');
if ~isempty(result_filename_override)
    cfg.result_filename = result_filename_override;
end

margin_values_override = getenv('PWQ_REV1_MARGIN_VALUES');
if ~isempty(margin_values_override)
    cfg.margin_values = parse_numeric_list(margin_values_override, 'PWQ_REV1_MARGIN_VALUES');
end

out_dir = fileparts(mfilename('fullpath'));
release_root = find_release_root(out_dir);
addpath(fullfile(release_root, 'data', 'double_integrator'));
configure_mosek_path(cfg);
result_csv = fullfile(ensure_generated_output_dir(release_root), cfg.result_filename);

reg = answer_regions();
scale_ctx = build_scale_context(reg, cfg);

fprintf('Running Rev1 fixed-Ts regional Lyapunov SDP scan.\n');
fprintf('Multiplier modes: %s, margins = %s, solver = %s\n', ...
    strjoin(cfg.multiplier_modes, ', '), mat2str(cfg.margin_values), cfg.solver);
fprintf('State scale: [%g, %g]\n', cfg.state_scale(1), cfg.state_scale(2));
fprintf('Row normalization: %s, redundancy: %s\n', cfg.row_normalization_mode, cfg.redundancy_mode);
fprintf('Origin regions: [%s]\n', strjoin(string(scale_ctx.origin_regions), ' '));
fprintf('Method-1 forward invariance status: %s\n', scale_ctx.forward_invariance_status);
fprintf('Partition complement max margin: %.3g over %d disjunctive cases\n', ...
    scale_ctx.partition_proof.max_complement_margin, scale_ctx.partition_proof.checked_cases);

rows = {};
total_runs = numel(cfg.margin_values) * numel(cfg.multiplier_modes) * numel(cfg.sample_times);
run_idx = 0;
for ee = 1:numel(cfg.margin_values)
    for mm = 1:numel(cfg.multiplier_modes)
        cfg_run = cfg;
        cfg_run.epsilon = cfg.margin_values(ee);
        cfg_run.alpha = cfg.margin_values(ee);
        cfg_run.multiplier_mode = cfg.multiplier_modes(mm);
        for kk = 1:numel(cfg.sample_times)
            Ts = cfg.sample_times(kk);
            run_idx = run_idx + 1;
            fprintf('\n[%d/%d] margin = %.3g, mode = %s, Ts = %.6g\n', ...
                run_idx, total_runs, cfg_run.epsilon, cfg_run.multiplier_mode, Ts);
            row = solve_one_Ts(Ts, reg, scale_ctx, cfg_run);
            rows{end+1} = row; %#ok<SAGROW>
            fprintf('  status: %s, transitions: %d, min_cheby_radius: %.3g, postcheck_max_eig: %.3g, beta: %.3g, solve_time: %.2fs\n', ...
                row.certificate_status, row.reachable_transition_count, ...
                row.min_transition_cheby_radius, row.postcheck_max_eig, ...
                row.beta_upper_bound, row.solve_time_sec);
        end
    end
end

T = struct2table([rows{:}]);
writetable(T, result_csv);
fprintf('\nWrote results:\n  %s\n', result_csv);


function row = solve_one_Ts(Ts, reg_law, scale_ctx, cfg)
n = 2;
M = numel(reg_law);
E = blkdiag(eye(n), 0);
reg_poly = scale_ctx.reg_poly;

[Phi, q, Abar] = build_sample_maps(reg_law, Ts, cfg.state_scale);
[partition_miss_count, partition_multi_count] = partition_grid_diagnostic(reg_poly, cfg.grid_x1, cfg.grid_x2);

[transitions, min_radius, degenerate_count, near_degenerate_count, ...
    transition_zero_rows, transition_redundant_rows, ambiguous_count] = ...
    enumerate_transitions(reg_poly, Phi, q, cfg);

yalmip('clear');
F = [];
postcheck_exprs = {};
L = cell(1, M);
P = cell(1, M);

for i = 1:M
    P{i} = sdpvar(n, n, 'symmetric');
    if ismember(i, scale_ctx.origin_regions)
        r = zeros(n, 1);
        c = 0;
    else
        r = sdpvar(n, 1);
        c = sdpvar(1, 1);
    end
    L{i} = [P{i}, r; r', c];
end

for i = 1:M
    Qpos = L{i} - cfg.epsilon * E;
    [F, postcheck_expr] = add_poly_quadratic_nonnegative_constraint(F, Qpos, ...
        reg_poly(i).H, reg_poly(i).b, cfg);
    postcheck_exprs{end+1} = postcheck_expr; %#ok<AGROW>
end

for k = 1:numel(transitions)
    i = transitions(k).i;
    j = transitions(k).j;
    Gij = transitions(k).G;
    wij = transitions(k).w;
    Qdec = L{i} - Abar{i}' * L{j} * Abar{i} - cfg.alpha * E;
    [F, postcheck_expr] = add_poly_quadratic_nonnegative_constraint(F, Qdec, Gij, wij, cfg);
    postcheck_exprs{end+1} = postcheck_expr; %#ok<AGROW>
end

ops = make_solver_options(cfg);
tic;
sol = optimize(F, [], ops);
solve_time = toc;

sdp_feasible = (sol.problem == 0);
[postcheck_max_eig, postcheck_pass] = postcheck_lmis(postcheck_exprs, cfg.postcheck_tol, sdp_feasible);

dense = empty_dense_postcheck();
if cfg.dense_postcheck_enabled && Ts <= cfg.dense_postcheck_max_Ts + 10*eps && sdp_feasible && postcheck_pass
    t_dense = tic;
    dense = dense_grid_postcheck(L, reg_poly, Phi, q, transitions, cfg);
    dense.elapsed_time_sec = toc(t_dense);
end

[beta_upper_bound, beta_postcheck_max_eig, beta_postcheck_pass, ...
    beta_problem_code, beta_yalmip_info] = ...
    postcheck_beta_upper_bound(L, reg_poly, cfg, sdp_feasible && postcheck_pass);

if sdp_feasible && postcheck_pass && beta_postcheck_pass && scale_ctx.forward_invariance_verified
    certificate_status = "certified";
elseif sdp_feasible && postcheck_pass && beta_postcheck_pass
    certificate_status = "conditional_on_forward_invariance";
elseif sdp_feasible && postcheck_pass
    certificate_status = "not_certified_beta_postcheck_failed";
elseif sdp_feasible
    certificate_status = "not_certified_postcheck_failed";
else
    certificate_status = "not_certified";
end

row = struct();
row.Ts = Ts;
row.epsilon = cfg.epsilon;
row.alpha = cfg.alpha;
row.multiplier_mode = string(cfg.multiplier_mode);
row.solver = string(cfg.solver);
row.mosek_primal_feas_tol = cfg.mosek_primal_feas_tol;
row.mosek_dual_feas_tol = cfg.mosek_dual_feas_tol;
row.mosek_rel_gap_tol = cfg.mosek_rel_gap_tol;
row.mosek_max_iter = cfg.mosek_max_iter;
row.mosek_num_threads = cfg.mosek_num_threads;
row.state_scale_s1 = cfg.state_scale(1);
row.state_scale_s2 = cfg.state_scale(2);
row.row_normalization_mode = cfg.row_normalization_mode;
row.redundancy_mode = cfg.redundancy_mode;
row.equality_multiplier_mode = cfg.equality_multiplier_mode;
row.boundary_continuity_mode = cfg.boundary_continuity_mode;
row.preprocess_safety_status = cfg.preprocess_safety_status;
row.sdp_problem_code = sol.problem;
row.sdp_yalmip_info = string(sol.info);
row.sdp_feasible = sdp_feasible;
row.postcheck_max_eig = postcheck_max_eig;
row.postcheck_pass = postcheck_pass;
row.postcheck_tol = cfg.postcheck_tol;
row.dense_grid_coordinate_system = string(cfg.dense_grid_coordinate_system);
row.dense_grid_scaled_x1_min = min(cfg.dense_grid_scaled_x1);
row.dense_grid_scaled_x1_max = max(cfg.dense_grid_scaled_x1);
row.dense_grid_scaled_x1_count = numel(cfg.dense_grid_scaled_x1);
row.dense_grid_scaled_x2_min = min(cfg.dense_grid_scaled_x2);
row.dense_grid_scaled_x2_max = max(cfg.dense_grid_scaled_x2);
row.dense_grid_scaled_x2_count = numel(cfg.dense_grid_scaled_x2);
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
row.beta_upper_bound = beta_upper_bound;
row.beta_postcheck_max_eig = beta_postcheck_max_eig;
row.beta_postcheck_pass = beta_postcheck_pass;
row.beta_sdp_problem_code = beta_problem_code;
row.beta_sdp_yalmip_info = beta_yalmip_info;
row.certificate_status = certificate_status;
row.forward_invariance_status = scale_ctx.forward_invariance_status;
row.full_space_partition_verified = scale_ctx.forward_invariance_verified;
row.partition_complement_max_margin = scale_ctx.partition_proof.max_complement_margin;
row.partition_complement_checked_cases = scale_ctx.partition_proof.checked_cases;
row.partition_complement_witness_x1 = scale_ctx.partition_proof.best_x(1);
row.partition_complement_witness_x2 = scale_ctx.partition_proof.best_x(2);
row.partition_complement_best_choice = strjoin(string(scale_ctx.partition_proof.best_choice), ' ');
row.partition_grid_miss_count = partition_miss_count;
row.partition_grid_multi_count = partition_multi_count;
row.origin_region_count = numel(scale_ctx.origin_regions);
row.origin_regions = strjoin(string(scale_ctx.origin_regions), ' ');
row.reachable_transition_count = numel(transitions);
row.ambiguous_transition_count = ambiguous_count;
row.degenerate_transition_count = degenerate_count;
row.near_degenerate_transition_count = near_degenerate_count;
row.min_transition_cheby_radius = min_radius;
row.zero_row_removed_count = scale_ctx.region_zero_row_removed_count + transition_zero_rows;
row.removed_region_constraint_count = scale_ctx.removed_region_constraint_count;
row.removed_transition_constraint_count = transition_redundant_rows;
row.solve_time_sec = solve_time;
end


function [beta_upper_bound, beta_postcheck_max_eig, beta_postcheck_pass, ...
    beta_problem_code, beta_yalmip_info] = ...
    postcheck_beta_upper_bound(L, reg_poly, cfg, enabled)
% Posterior check for the upper comparison bound in paper Eq. (6a).
% Numeric L_i are fixed, so this auxiliary SDP does not change the main search.
if ~enabled
    beta_upper_bound = NaN;
    beta_postcheck_max_eig = NaN;
    beta_postcheck_pass = false;
    beta_problem_code = NaN;
    beta_yalmip_info = "not_run";
    return;
end

M = numel(reg_poly);
n = size(L{1}, 1) - 1;
E = blkdiag(eye(n), 0);
L_value = cell(1, M);
for i = 1:M
    Li = value(L{i});
    Li = (Li + Li') / 2;
    if any(~isfinite(Li), 'all')
        beta_upper_bound = NaN;
        beta_postcheck_max_eig = inf;
        beta_postcheck_pass = false;
        beta_problem_code = NaN;
        beta_yalmip_info = "nonfinite_L_value";
        return;
    end
    L_value{i} = Li;
end

yalmip('clear');
beta = sdpvar(1, 1);
F = [beta >= cfg.epsilon];
postcheck_exprs = {};
for i = 1:M
    Qupper = beta * E - L_value{i};
    [F, postcheck_expr] = add_poly_quadratic_nonnegative_constraint(F, Qupper, ...
        reg_poly(i).H, reg_poly(i).b, cfg);
    postcheck_exprs{end+1} = postcheck_expr; %#ok<AGROW>
end

ops = make_solver_options(cfg);
sol = optimize(F, beta, ops);
beta_problem_code = sol.problem;
beta_yalmip_info = string(sol.info);

if sol.problem == 0
    beta_upper_bound = value(beta);
else
    beta_upper_bound = NaN;
end

[beta_postcheck_max_eig, lmi_pass] = postcheck_lmis(postcheck_exprs, cfg.postcheck_tol, sol.problem == 0);
beta_postcheck_pass = lmi_pass && isfinite(beta_upper_bound) && beta_upper_bound > 0;
end


function [F, postcheck_expr] = add_poly_quadratic_nonnegative_constraint(F, Q, G, w, cfg)
% Adds sufficient constraints for z' Q z >= 0 over {x | G x <= w},
% where z = [x; 1].
nz = size(Q, 1);
m = size(G, 1);

switch string(cfg.multiplier_mode)
    case "scalar"
        lambda = sdpvar(m, 1);
        Mmat = Q;
        for ell = 1:m
            g = G(ell, :)';
            Sell = [zeros(nz-1, nz-1), -g/2; -g'/2, w(ell)];
            Mmat = Mmat - lambda(ell) * Sell;
        end
        F = [F, lambda >= 0, Mmat >= 0]; %#ok<AGROW>
        postcheck_expr = -Mmat;

    case "matrix"
        Gamma = [-G, w; zeros(1, nz-1), 1];
        Z = sdpvar(m+1, m+1, 'symmetric');
        Mmat = Q - Gamma' * Z * Gamma;
        F = [F, Z(:) >= 0, Mmat >= 0]; %#ok<AGROW>
        postcheck_expr = -Mmat;

    otherwise
        error('Unknown multiplier mode: %s', cfg.multiplier_mode);
end
end


function [max_eig, pass] = postcheck_lmis(exprs, tol, enabled)
% Checks LMI matrices after solving. Each expression is expected to satisfy
% expr <= 0, so max_eig <= tol is the acceptance criterion.
if ~enabled
    max_eig = NaN;
    pass = false;
    return;
end

max_eig = -inf;
for k = 1:numel(exprs)
    M = value(exprs{k});
    M = (M + M') / 2;
    if any(~isfinite(M), 'all')
        max_eig = inf;
        pass = false;
        return;
    end
    max_eig = max(max_eig, max(eig(M)));
end

if isempty(exprs)
    max_eig = NaN;
    pass = false;
else
    pass = max_eig <= tol;
end
end


function dense = dense_grid_postcheck(L, reg_poly, Phi, q, transitions, cfg)
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
        xh = [xv; 1];
        point_count = point_count + 1;

        hit_regions = [];
        for i = 1:M
            if all(reg_poly(i).H * xv <= reg_poly(i).b + cfg.lp_tol)
                hit_regions(end+1) = i; %#ok<AGROW>
                Li = value(L{i});
                val = xh' * Li * xh - cfg.epsilon * (xv' * xv);
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
                xph = [xp; 1];
                Li = value(L{i});
                Lj = value(L{j});
                val = xh' * Li * xh - xph' * Lj * xph - cfg.alpha * (xv' * xv);
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


function values = parse_numeric_list(raw_value, name)
cleaned = regexprep(strtrim(raw_value), '[,;]+', ' ');
values = sscanf(cleaned, '%f')';
if isempty(values) || any(~isfinite(values))
    error('Invalid %s: %s', name, raw_value);
end
end


function status = build_preprocess_safety_status(cfg)
parts = strings(1, 0);
if cfg.row_normalization_mode == "unit_normal"
    parts(end+1) = "equivalent_unit_row_normalization"; %#ok<AGROW>
end
if cfg.redundancy_mode == "remove_redundant"
    parts(end+1) = "outer_or_equal_redundancy_removal"; %#ok<AGROW>
end
if isempty(parts)
    status = "raw_no_preprocessing";
else
    status = strjoin(parts, "+");
end
end


function scale_ctx = build_scale_context(reg, cfg)
T = diag(cfg.state_scale);
scale_ctx = struct();
scale_ctx.state_scale_matrix = T;
[scale_ctx.reg_poly, scale_ctx.region_zero_row_removed_count, ...
    scale_ctx.removed_region_constraint_count] = ...
    scale_region_polytopes(reg, T, cfg);
[~, scale_ctx.origin_regions] = find_origin_regions(scale_ctx.reg_poly);
scale_ctx.partition_proof = verify_full_state_space_partition(scale_ctx.reg_poly, cfg.lp_tol);
if scale_ctx.partition_proof.full_space_covered
    scale_ctx.forward_invariance_verified = true;
    scale_ctx.forward_invariance_status = "verified_full_state_space_partition";
else
    scale_ctx.forward_invariance_verified = false;
    scale_ctx.forward_invariance_status = "not_verified";
end
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


function [Phi, q, Abar] = build_sample_maps(reg, Ts, state_scale)
M = numel(reg);
Phi = cell(1, M);
q = cell(1, M);
Abar = cell(1, M);
A0 = [1 Ts; 0 1];
B0 = [Ts^2/2; Ts];
T = diag(state_scale);

for i = 1:M
    Phi_original = A0 + B0 * reg(i).K;
    q_original = B0 * reg(i).k;
    Phi{i} = T \ (Phi_original * T);
    q{i} = T \ q_original;
    Abar{i} = [Phi{i}, q{i}; 0 0 1];
end
end


function [Gout, wout, zero_removed, removed_redundant, is_empty] = preprocess_polyhedron(G, w, cfg)
% Keeps the represented polyhedron unchanged except for LP-certified
% redundant rows. If numerical tolerance has an effect, row removal is an
% outer/equal approximation direction, never an inner approximation.
row_norms = vecnorm(G, 2, 2);
keep = row_norms >= cfg.zero_row_tol;
zero_rows = ~keep;
zero_removed = 0;
removed_redundant = 0;
is_empty = false;

if any(zero_rows)
    if any(w(zero_rows) < -cfg.zero_row_tol)
        Gout = G;
        wout = w;
        is_empty = true;
        return;
    end
    zero_removed = nnz(zero_rows);
end

Gout = G(keep, :);
wout = w(keep);
row_norms = row_norms(keep);

switch string(cfg.row_normalization_mode)
    case "none"

    case "unit_normal"
        if ~isempty(Gout)
            Gout = Gout ./ row_norms;
            wout = wout ./ row_norms;
        end

    otherwise
        error('Unknown row normalization mode: %s', cfg.row_normalization_mode);
end

if cfg.redundancy_mode == "remove_redundant"
    [Gout, wout, removed_redundant] = remove_redundant_constraints(Gout, wout, cfg.redundancy_tol);
elseif cfg.redundancy_mode ~= "keep_redundant"
    error('Unknown redundancy mode: %s', cfg.redundancy_mode);
end
end


function [Gkeep, wkeep, removed_count] = remove_redundant_constraints(G, w, tol)
m = size(G, 1);
if m <= 1
    Gkeep = G;
    wkeep = w;
    removed_count = 0;
    return;
end

keep = true(m, 1);
opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');

for ell = 1:m
    other = true(m, 1);
    other(ell) = false;
    f = -G(ell, :)';
    [~, fval, exitflag] = linprog(f, G(other, :), w(other), [], [], [], [], opts);
    if exitflag == 1
        max_violation = -fval - w(ell);
        if max_violation <= tol
            keep(ell) = false;
        end
    end
end

Gkeep = G(keep, :);
wkeep = w(keep);
removed_count = nnz(~keep);
end


function [transitions, min_radius, degenerate_count, near_degenerate_count, ...
    zero_removed_total, removed_redundant_total, ambiguous_count] = enumerate_transitions(reg, Phi, q, cfg)
M = numel(reg);
transitions = struct('i', {}, 'j', {}, 'G', {}, 'w', {}, 'cheby_radius', {});
radii = [];
degenerate_count = 0;
near_degenerate_count = 0;
zero_removed_total = 0;
removed_redundant_total = 0;
ambiguous_count = 0;

for i = 1:M
    for j = 1:M
        Gij = [reg(i).H; reg(j).H * Phi{i}];
        wij = [reg(i).b; reg(j).b - reg(j).H * q{i}];
        [Gij, wij, zero_removed, removed_redundant, is_empty] = preprocess_polyhedron(Gij, wij, cfg);
        zero_removed_total = zero_removed_total + zero_removed;
        removed_redundant_total = removed_redundant_total + removed_redundant;
        if is_empty
            continue;
        end
        [is_feas, radius, is_ambiguous] = poly_feasibility_and_cheby_radius(Gij, wij, cfg.lp_tol);
        if is_feas
            if is_ambiguous
                ambiguous_count = ambiguous_count + 1;
            end
            transitions(end+1).i = i; %#ok<AGROW>
            transitions(end).j = j;
            transitions(end).G = Gij;
            transitions(end).w = wij;
            transitions(end).cheby_radius = radius;
            radii(end+1) = radius; %#ok<AGROW>
            if radius <= max(cfg.lp_tol, 1e-7)
                degenerate_count = degenerate_count + 1;
            end
            if radius <= cfg.near_degenerate_radius_tol
                near_degenerate_count = near_degenerate_count + 1;
            end
        end
    end
end

if isempty(radii)
    min_radius = NaN;
else
    min_radius = min(radii);
end
end


function [is_feas, radius, is_ambiguous] = poly_feasibility_and_cheby_radius(G, w, tol)
n = size(G, 2);
is_ambiguous = false;
f = [zeros(n,1); -1];
norms = vecnorm(G, 2, 2);
A = [G, norms];
b = w;
lb = [-inf(n,1); 0];

opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
[~, fval, exitflag] = linprog(f, A, b, [], [], lb, [], opts);

is_feas = exitflag == 1 || exitflag == -3;
if is_feas
    if exitflag == -3
        radius = inf;
    else
        radius = max(0, -fval);
    end
else
    radius = NaN;
end

if ~is_feas
    % Conservative fallback: if numerical pure feasibility cannot rule out
    % the transition, keep it and mark it ambiguous.
    [~, ~, exitflag2] = linprog(zeros(n,1), G, w + tol, [], [], [], [], opts);
    is_feas = exitflag2 == 1 || exitflag2 == -3;
    if is_feas
        radius = 0;
        is_ambiguous = true;
    end
end
end


function [miss_count, multi_count] = partition_grid_diagnostic(reg, x1_grid, x2_grid)
miss_count = 0;
multi_count = 0;
for a = 1:numel(x1_grid)
    for b = 1:numel(x2_grid)
        x = [x1_grid(a); x2_grid(b)];
        hits = 0;
        for i = 1:numel(reg)
            if all(reg(i).H * x <= reg(i).b + 1e-9)
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


function configure_mosek_path(cfg)
if ~isfield(cfg, 'mosek_path') || isempty(cfg.mosek_path)
    return;
end
mosek_path = char(cfg.mosek_path);
if isfolder(mosek_path) && ~any(strcmpi(strsplit(path, pathsep), mosek_path))
    addpath(mosek_path);
end
end


function root = find_release_root(start_dir)
root = start_dir;
while true
    if exist(fullfile(root, 'data', 'double_integrator', 'answer_regions.m'), 'file') == 2
        return;
    end
    parent = fileparts(root);
    if strcmp(parent, root)
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


function ops = make_solver_options(cfg)
ops = sdpsettings('solver', char(cfg.solver), 'verbose', cfg.verbose);
ops.cachesolvers = 1;

if strcmpi(string(cfg.solver), "mosek")
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = cfg.mosek_primal_feas_tol;
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = cfg.mosek_dual_feas_tol;
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = cfg.mosek_rel_gap_tol;
    ops.mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS = cfg.mosek_max_iter;
    ops.mosek.MSK_IPAR_NUM_THREADS = cfg.mosek_num_threads;
end
end


function [idx, hits] = find_origin_regions(reg)
x0 = [0; 0];
hits = [];
for i = 1:numel(reg)
    if all(reg(i).H * x0 <= reg(i).b + 1e-9)
        hits(end+1) = i; %#ok<AGROW>
    end
end

if isempty(hits)
    error('Origin is not covered by any explicit-MPC region.');
end
idx = hits(1);
end


function proof = verify_full_state_space_partition(reg, tol)
% Verifies whether union_i {x | H_i x <= b_i} covers R^2.
M = numel(reg);
row_counts = arrayfun(@(r) size(r.H, 1), reg);
choices = ones(1, M);
checked_cases = 0;
max_margin = -inf;
best_choice = choices;
best_x = [NaN; NaN];

opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');

while true
    checked_cases = checked_cases + 1;
    A = zeros(M, 3);
    b = zeros(M, 1);
    for i = 1:M
        row = choices(i);
        h = reg(i).H(row, :);
        bi = reg(i).b(row);
        A(i, :) = [-h, 1];
        b(i) = -bi;
    end

    f = [0; 0; -1];
    [sol, fval, exitflag] = linprog(f, A, b, [], [], [], [], opts);

    if exitflag == -3
        margin = inf;
    elseif exitflag == 1
        margin = -fval;
    else
        margin = -inf;
    end

    if margin > max_margin
        max_margin = margin;
        best_choice = choices;
        if exist('sol', 'var') && numel(sol) >= 2
            best_x = sol(1:2);
        end
    end

    if margin > tol
        break;
    end

    pos = M;
    while pos >= 1
        choices(pos) = choices(pos) + 1;
        if choices(pos) <= row_counts(pos)
            break;
        end
        choices(pos) = 1;
        pos = pos - 1;
    end
    if pos == 0
        break;
    end
end

proof = struct();
proof.full_space_covered = max_margin <= tol;
proof.max_complement_margin = max_margin;
proof.checked_cases = checked_cases;
proof.best_choice = best_choice;
proof.best_x = best_x;
end
