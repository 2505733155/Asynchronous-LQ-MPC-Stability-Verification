%% KKT-lifted SOS certificate for sampled standard MPC
% First implementation of the Korda-Jones KKT/SOS route for the online QP
% standard MPC defined by the package-local seven-region PWA data.
%
% By default the script searches local compact-domain sampled-point SOS
% decrease certificates. When KKT_SOS_GLOBAL_ATTEMPT=1, it runs a direct
% noncompact full-space comparison attempt. Neither mode proves continuous
% inter-sample stability.

clear; clc;

cfg = make_cfg();
script_dir = get_script_dir();
release_root = find_release_root(script_dir);
addpath(script_dir);
addpath(fullfile(release_root, 'data', 'double_integrator'));
ensure_yalmip_export_precedence();
configure_mosek_path(cfg);
ensure_matlab_linprog_precedence();
output_dir = ensure_generated_output_dir(release_root);

if cfg.global_attempt_mode
    out_csv = fullfile(output_dir, ['kkt_lifted_sos_standard_mpc_global_attempt_results' char(cfg.output_suffix) '.csv']);
    out_mat = fullfile(output_dir, ['kkt_lifted_sos_standard_mpc_global_attempt_certificate' char(cfg.output_suffix) '.mat']);
    out_coeff_debug_csv = fullfile(output_dir, ['kkt_lifted_sos_global_attempt_coeff_debug' char(cfg.output_suffix) '.csv']);
else
    out_csv = fullfile(output_dir, ['kkt_lifted_sos_standard_mpc_results' char(cfg.output_suffix) '.csv']);
    out_mat = fullfile(output_dir, ['kkt_lifted_sos_standard_mpc_certificate' char(cfg.output_suffix) '.mat']);
    out_coeff_debug_csv = fullfile(output_dir, ['kkt_lifted_sos_coeff_debug' char(cfg.output_suffix) '.csv']);
end

fprintf('KKT-lifted SOS standard MPC certificate generator\n');
fprintf('Output CSV: %s\n', out_csv);
fprintf('Optional MAT output: %s\n', out_mat);
if cfg.coeff_debug_enabled
    fprintf('Coefficient debug CSV: %s\n', out_coeff_debug_csv);
end
fprintf('Run label: %s\n', cfg.run_label);

mpc = build_standard_mpc_data();
precheck = precheck_standard_mpc_data(mpc, cfg);
cases = build_case_list(cfg);
rows = [];
coeff_debug_rows = [];
cert_bundle = struct();

for kk = 1:numel(cases)
    this_case = cases(kk);
    fprintf('\n[KKT/SOS] Ts=%.12g, d_V=%d, eps=%.3g, alpha=%.3g\n', ...
        this_case.Ts, this_case.dV, cfg.epsilon, cfg.alpha);
    [row, cert] = solve_one_case(this_case.Ts, this_case.dV, mpc, precheck, cfg);
    rows = append_row(rows, row);
    coeff_debug_rows = append_coeff_debug_rows(coeff_debug_rows, cert, row, cfg);
    cert_bundle.(this_case.key) = cert;
    print_case_summary(row);
end

T = struct2table(rows);
writetable(T, out_csv);
if cfg.coeff_debug_enabled && ~isempty(coeff_debug_rows)
    writetable(struct2table(coeff_debug_rows), out_coeff_debug_csv);
end
if strcmp(getenv('AMPC_SAVE_MAT'), '1')
    save(out_mat, 'cert_bundle', 'cfg', 'mpc', 'precheck');
end

disp(T(:, {'Ts','d_V','certificate_status','solver_problem_code', ...
    'sos_constraint_count','gram_min_eig','coeff_match_max_abs', ...
    'qp_collision_min_decay','roa_gamma','solve_time_sec'}));


function cfg = make_cfg()
cfg.run_label = "local_box_Ts_1_0p1_dV2_4_exact_kkt";
cfg.sample_times = [1.0, 0.1];
cfg.dV_values = [2, 4];
cfg.use_parallel = false;

cfg.epsilon = 1e-6;
cfg.alpha = 1e-6;
cfg.domain_x1 = [-15, 15];
cfg.domain_x2 = [-6, 6];
cfg.domain_profile = "baseline_box";
cfg.output_suffix = "";
cfg.kkt_mode = "exact";
cfg.relaxed_delta = NaN;
cfg.certificate_claim_mode = "sos_decrease_no_invariance";
cfg.require_invariance_for_status = false;
cfg.global_attempt_mode = false;
cfg.noncompact_lifted_set = false;
cfg.analysis_domain = "local_box";

cfg.solver = "mosek";
cfg.mosek_path = getenv("AMPC_MOSEK_PATH");
cfg.verbose = 0;
cfg.yalmip_sos_model = 0;
cfg.yalmip_sos_scale = 1;
cfg.yalmip_sos_newton = 0;
cfg.stationarity_normalization = "none";
cfg.lift_variable_scaling = "none";
cfg.fixed_state_scale = [20; 5];
cfg.coeff_debug_enabled = false;
cfg.coeff_debug_random_points = 20;
cfg.mosek_primal_feas_tol = 1e-8;
cfg.mosek_dual_feas_tol = 1e-8;
cfg.mosek_rel_gap_tol = 1e-8;
cfg.mosek_max_iter = 200;
cfg.mosek_num_threads = 1;

cfg.coefficient_bound = 1e3;
cfg.sos_multiplier_degree = 2;
cfg.domain_multiplier_degree = 2;
cfg.stationarity_multiplier_degree_dV2 = 1;
cfg.stationarity_multiplier_degree_dV4 = 3;
cfg.complementarity_multiplier_degree = 2;
cfg.max_total_gram_dim = 2500;

cfg.gram_eig_tol = 1e-7;
cfg.coeff_match_tol = 1e-2;
cfg.sos_oracle_tol_factor = 10;
cfg.yalmip_check_tol = 1e-6;
cfg.qp_collision_tol = 1e-7;
cfg.qp_collision_tol_factor = 10;
cfg.roa_bisection_steps = 8;
cfg.roa_gamma_floor = 1e-8;

cfg.collision_grid_x1 = linspace(-15, 15, 9);
cfg.collision_grid_x2 = linspace(-6, 6, 7);
cfg.collision_random_points = 100;
cfg.collision_random_seed = 7;
cfg.origin_probe_radius = 1e-5;
cfg.global_sanity_x1 = [-50, 50];
cfg.global_sanity_x2 = [-20, 20];

if strcmp(getenv('KKT_SOS_SMOKE'), '1')
    cfg.run_label = "smoke_Ts_1_dV2_exact_kkt";
    cfg.sample_times = 1.0;
    cfg.dV_values = 2;
    cfg.collision_random_points = 20;
    cfg.roa_bisection_steps = 4;
end

if strcmp(getenv('KKT_SOS_COEFF_DEBUG'), '1')
    cfg.coeff_debug_enabled = true;
    cfg.run_label = "coeff_debug_Ts_0p1_1_dV2_exact_kkt";
    cfg.sample_times = [0.1, 1.0];
    cfg.dV_values = 2;
    cfg.collision_random_points = 20;
    cfg.roa_bisection_steps = 4;
end

if strcmp(getenv('KKT_SOS_GLOBAL_ATTEMPT'), '1')
    cfg.global_attempt_mode = true;
    cfg.noncompact_lifted_set = true;
    cfg.analysis_domain = "R2";
    cfg.run_label = "global_attempt_Ts_1_0p1_dV2_4_exact_kkt";
    cfg.sample_times = [1.0, 0.1];
    cfg.dV_values = [2, 4];
    cfg.certificate_claim_mode = "global_sampled_asymptotic";
    cfg.require_invariance_for_status = false;
    cfg.max_total_gram_dim = Inf;
    cfg.domain_profile = "global_attempt_noncompact";
    cfg.domain_x1 = cfg.global_sanity_x1;
    cfg.domain_x2 = cfg.global_sanity_x2;
    cfg.roa_bisection_steps = 0;
end

env_domain_profile = getenv('KKT_SOS_DOMAIN_PROFILE');
if ~isempty(env_domain_profile) && ~cfg.global_attempt_mode
    cfg.domain_profile = string(env_domain_profile);
    cfg = apply_domain_profile(cfg);
end

env_output_suffix = getenv('KKT_SOS_OUTPUT_SUFFIX');
if ~isempty(env_output_suffix)
    cfg.output_suffix = string(env_output_suffix);
end

env_sos_scale = getenv('KKT_SOS_SCALE');
if ~isempty(env_sos_scale)
    cfg.yalmip_sos_scale = str2double(env_sos_scale);
end
env_sos_model = getenv('KKT_SOS_MODEL');
if ~isempty(env_sos_model)
    cfg.yalmip_sos_model = str2double(env_sos_model);
end
env_stationarity = getenv('KKT_SOS_STATIONARITY_NORMALIZATION');
if ~isempty(env_stationarity)
    cfg.stationarity_normalization = string(env_stationarity);
end
env_lift_scaling = getenv('KKT_SOS_LIFT_SCALING');
if ~isempty(env_lift_scaling)
    cfg.lift_variable_scaling = string(env_lift_scaling);
end
if cfg.global_attempt_mode
    if strlength(strtrim(env_lift_scaling)) == 0
        cfg.lift_variable_scaling = "none";
    end
end

if cfg.coeff_debug_enabled
    cfg.run_label = cfg.run_label + "_scale" + string(cfg.yalmip_sos_scale) + ...
        "_model" + string(cfg.yalmip_sos_model) + "_stat_" + ...
        cfg.stationarity_normalization + "_lift_" + cfg.lift_variable_scaling + ...
        "_domain_" + cfg.domain_profile;
end
cfg = sync_collision_grid_to_domain(cfg);
end


function cfg = apply_domain_profile(cfg)
switch string(cfg.domain_profile)
    case "baseline_box"
        cfg.domain_x1 = [-15, 15];
        cfg.domain_x2 = [-6, 6];
    case "small_box"
        cfg.domain_x1 = [-5, 5];
        cfg.domain_x2 = [-3, 3];
        cfg.run_label = cfg.run_label + "_small_box_x1_5_x2_3";
    case "tiny_box"
        cfg.domain_x1 = [-3, 3];
        cfg.domain_x2 = [-2, 2];
        cfg.run_label = cfg.run_label + "_tiny_box_x1_3_x2_2";
    otherwise
        error('Unknown KKT_SOS_DOMAIN_PROFILE: %s', cfg.domain_profile);
end
end


function cfg = sync_collision_grid_to_domain(cfg)
cfg.collision_grid_x1 = linspace(cfg.domain_x1(1), cfg.domain_x1(2), 9);
cfg.collision_grid_x2 = linspace(cfg.domain_x2(1), cfg.domain_x2(2), 7);
end


function cases = build_case_list(cfg)
cases = struct('Ts', {}, 'dV', {}, 'key', {});
kk = 0;
for tt = 1:numel(cfg.sample_times)
    for dd = 1:numel(cfg.dV_values)
        kk = kk + 1;
        Ts = cfg.sample_times(tt);
        dV = cfg.dV_values(dd);
        cases(kk).Ts = Ts; %#ok<AGROW>
        cases(kk).dV = dV; %#ok<AGROW>
        cases(kk).key = matlab.lang.makeValidName(sprintf('Ts_%s_dV_%d', num2str(Ts), dV)); %#ok<AGROW>
    end
end
end


function [row, cert] = solve_one_case(Ts, dV, mpc, precheck, cfg)
t_case = tic;
row = base_row(Ts, dV, cfg, precheck);
cert = empty_cert(Ts, dV, cfg);

if ~precheck.tools_available
    row.certificate_status = status_solver_failed(cfg);
    row.failure_stage = "tool_precheck";
    row.solver_info = precheck.tool_info;
    row.solve_time_sec = toc(t_case);
    cert.status = row.certificate_status;
    return;
end

try
    plant = sampled_execution_model(Ts);
    bounds = build_lift_bounds(plant, mpc, cfg);
    row.R_lambda = bounds.R_lambda;
    row.R_eta = bounds.R_eta;
    row.execution_prediction_mismatch_norm = norm(plant.As - mpc.Ad, 'fro') + norm(plant.Bs - mpc.Bd, 'fro');
    row.lifted_variable_count = 14;

    scales = build_variable_scaling(bounds, cfg);
    [F, params, poly_meta] = build_kkt_sos_problem(Ts, dV, plant, mpc, bounds, scales, cfg);
    row.total_gram_dim_estimate = sum([poly_meta.sos_meta.gram_dim]);
    row.sos_constraint_count = nnz([poly_meta.sos_meta.is_sos_constraint]);

    if row.total_gram_dim_estimate > cfg.max_total_gram_dim
        row.certificate_status = "not_attempted_due_to_size";
        row.failure_stage = "size_guard";
        row.solver_info = sprintf('total_gram_dim_estimate=%g exceeds max_total_gram_dim=%g', ...
            row.total_gram_dim_estimate, cfg.max_total_gram_dim);
        row.solve_time_sec = toc(t_case);
        cert.status = row.certificate_status;
        cert.poly_meta = strip_runtime_handles(poly_meta);
        return;
    end

    ops = make_solver_options(cfg);
    [sol, sos_monomials, sos_grams, sos_residuals] = solvesos( ...
        F, [], ops, params, poly_meta.sos_candidate_monomials); %#ok<ASGLU>
    row.solver_problem_code = sol.problem;
    row.solver_info = string(sol.info);
    row.solve_time_sec = toc(t_case);

    if sol.problem ~= 0
        if sol.problem == 1
            row.certificate_status = status_inconclusive(cfg);
            row.failure_stage = "sos_infeasible";
        else
            row.certificate_status = status_solver_failed(cfg);
            row.failure_stage = "sos_solver";
        end
        cert.status = row.certificate_status;
        cert.poly_meta = strip_runtime_handles(poly_meta);
        return;
    end

    check_vals = check(F);
    row.yalmip_check_min = min(check_vals);
    row.yalmip_check_strict_pass = row.yalmip_check_min >= -cfg.yalmip_check_tol;
    row.yalmip_check_tol_pass = row.yalmip_check_min >= -cfg.sos_oracle_tol_factor * cfg.yalmip_check_tol;

    sos_oracle = independent_sos_postcheck_kkt(poly_meta, cfg, ...
        numel(sos_grams), sos_monomials, sos_grams, sos_residuals);
    row.gram_min_eig = sos_oracle.min_gram_eig;
    row.gram_worst_block = sos_oracle.worst_gram_block;
    row.coeff_match_max_abs = sos_oracle.max_abs_coeff_residual;
    row.coeff_match_worst_block = sos_oracle.worst_coeff_block;
    row.sos_postcheck_gram_handle_missing_count = sos_oracle.gram_handle_missing_count;
    row.sos_postcheck_gram_dim_mismatch_count = sos_oracle.gram_dim_mismatch_count;
    row.sos_oracle_strict_pass = sos_oracle.strict_pass;
    row.sos_oracle_tol_pass = sos_oracle.tol_pass;

    V_fixed = fixed_polynomial(poly_meta.V_coeff, poly_meta.V_monomials);
    collision = qp_collision_postcheck(V_fixed, poly_meta.x, plant, mpc, cfg, scales);
    row.qp_collision_point_count = collision.point_count;
    row.qp_collision_min_decay = collision.min_decay;
    row.qp_collision_strict_pass = collision.strict_pass;
    row.qp_collision_tol_pass = collision.tol_pass;
    row.qp_collision_xp_out_box_count = collision.xp_out_box_count;
    row.qp_collision_worst_x1 = collision.worst_x(1);
    row.qp_collision_worst_x2 = collision.worst_x(2);
    row.qp_collision_worst_u1 = collision.worst_u1;
    row.qp_collision_worst_xp1 = collision.worst_xp(1);
    row.qp_collision_worst_xp2 = collision.worst_xp(2);
    row.qp_collision_worst_Vx = collision.worst_Vx;
    row.qp_collision_worst_Vxp = collision.worst_Vxp;
    row.qp_collision_worst_xp_in_box = collision.worst_xp_in_box;

    if cfg.global_attempt_mode
        roa = struct('gamma', NaN, 'status', "not_applicable_global_attempt", ...
            'strict_pass', false, 'tol_pass', false);
        collision_roa = struct('point_count', 0, 'failed_solves', 0, ...
            'min_decay', NaN, 'strict_pass', false, 'tol_pass', false);
    else
        roa = find_roa_gamma(V_fixed, poly_meta.x, cfg, scales);
        collision_roa = qp_collision_roa_postcheck(V_fixed, poly_meta.x, plant, mpc, cfg, scales, roa.gamma);
    end
    row.roa_gamma = roa.gamma;
    row.roa_invariance_status = roa.status;
    row.roa_inclusion_strict_pass = roa.strict_pass;
    row.roa_inclusion_tol_pass = roa.tol_pass;
    row.qp_collision_roa_point_count = collision_roa.point_count;
    row.qp_collision_roa_min_decay = collision_roa.min_decay;
    row.qp_collision_roa_strict_pass = collision_roa.strict_pass;
    row.qp_collision_roa_tol_pass = collision_roa.tol_pass;

    strict_pass = row.yalmip_check_strict_pass && row.sos_oracle_strict_pass;
    tol_pass = row.yalmip_check_tol_pass && row.sos_oracle_tol_pass;

    if cfg.global_attempt_mode
        global_strict_pass = strict_pass && row.qp_collision_strict_pass;
        global_tol_pass = tol_pass && row.qp_collision_tol_pass;
        if global_strict_pass
            row.certificate_status = "certified_global_sampled_asymptotic";
            row.failure_stage = "";
        elseif global_tol_pass
            row.certificate_status = "tol_certified_global_sampled_asymptotic";
            row.failure_stage = "tolerance_only";
        elseif strict_pass || tol_pass
            row.certificate_status = "inconclusive_global_attempt";
            row.failure_stage = "quadprog_sanity_failed";
        else
            row.certificate_status = "inconclusive_global_attempt";
            row.failure_stage = "sos_postcheck_failed";
        end
    elseif cfg.require_invariance_for_status
        strict_pass = strict_pass && row.qp_collision_strict_pass && row.roa_inclusion_strict_pass;
        tol_pass = tol_pass && row.qp_collision_tol_pass && row.roa_inclusion_tol_pass;
        if strict_pass
            row.certificate_status = "certified_asymptotic";
            row.failure_stage = "";
        elseif tol_pass
            row.certificate_status = "tol_certified_asymptotic";
            row.failure_stage = "tolerance_only";
        else
            row.certificate_status = "inconclusive";
            row.failure_stage = "postcheck_failed";
        end
    else
        if strict_pass
            row.certificate_status = "certified_sos_decrease_no_invariance";
            row.failure_stage = "";
        elseif tol_pass
            row.certificate_status = "tol_certified_sos_decrease_no_invariance";
            row.failure_stage = "tolerance_only";
        else
            row.certificate_status = "inconclusive";
            row.failure_stage = "sos_postcheck_failed";
        end
    end

    cert.status = row.certificate_status;
    cert.row = row;
    cert.V_coefficients = value(poly_meta.V_coeff);
    cert.bounds = bounds;
    cert.scales = scales;
    cert.sos_oracle = sos_oracle;
    cert.collision = collision;
    cert.collision_roa = collision_roa;
    cert.roa = roa;
    cert.poly_meta = strip_runtime_handles(poly_meta);
catch ME
    row.certificate_status = status_solver_failed(cfg);
    row.failure_stage = "exception";
    row.solver_info = string(getReport(ME, 'basic', 'hyperlinks', 'off'));
    row.solve_time_sec = toc(t_case);
    cert.status = row.certificate_status;
    cert.exception = row.solver_info;
end
end


function status = status_inconclusive(cfg)
if cfg.global_attempt_mode
    status = "inconclusive_global_attempt";
else
    status = "inconclusive";
end
end


function status = status_solver_failed(cfg)
if cfg.global_attempt_mode
    status = "solver_failed_global_attempt";
else
    status = "solver_failed";
end
end


function [F, params, poly_meta] = build_kkt_sos_problem(Ts, dV, plant, mpc, bounds, scales, cfg)
x_model = sdpvar(2, 1);
U_model = sdpvar(2, 1);
lub_model = sdpvar(2, 1);
llb_model = sdpvar(2, 1);
Un_model = sdpvar(2, 1);
eub_model = sdpvar(2, 1);
elb_model = sdpvar(2, 1);

x = scales.x(:) .* x_model;
U = scales.U(:) .* U_model;
lub = scales.lambda * lub_model;
llb = scales.lambda * llb_model;
Un = scales.U(:) .* Un_model;
eub = scales.eta * eub_model;
elb = scales.eta * elb_model;

lambda_model = [lub_model; llb_model];
eta_model = [eub_model; elb_model];
lambda = [lub; llb];
eta = [eub; elb];
z_current = [x_model; U_model; lambda_model];
z_lift = [x_model; U_model; lambda_model; Un_model; eta_model];
xplus = plant.As * x + plant.Bs * U(1);
xplus_model = xplus ./ scales.x(:);

[V, cV, vmon] = polynomial(x_model, dV);
params = cV;
F = [cV <= cfg.coefficient_bound, cV >= -cfg.coefficient_bound, ...
    replace(V, x_model, [0; 0]) == 0]; %#ok<NBRAK>
gradV = jacobian(V, x_model);
F = [F, replace(gradV(1), x_model, [0; 0]) == 0, ...
    replace(gradV(2), x_model, [0; 0]) == 0];

poly_meta = struct();
poly_meta.x = x_model;
poly_meta.x_physical = x;
poly_meta.U = U_model;
poly_meta.lambda = lambda_model;
poly_meta.Un = Un_model;
poly_meta.eta = eta_model;
poly_meta.z_current = z_current;
poly_meta.z_lift = z_lift;
poly_meta.V = V;
poly_meta.V_coeff = cV;
poly_meta.V_monomials = vmon;
poly_meta.dV = dV;
poly_meta.Ts = Ts;
poly_meta.scales = scales;
poly_meta.sos_meta = empty_sos_meta_kkt();
poly_meta.sos_candidate_monomials = {};

[h_cur, hdeg_cur, g_cur] = kkt_polynomials(x, U, lub, llb, ...
    lub_model, llb_model, mpc, cfg, scales.lambda);
[h_next, hdeg_next, g_next] = kkt_polynomials(xplus, Un, eub, elb, ...
    eub_model, elb_model, mpc, cfg, scales.eta);
if cfg.global_attempt_mode
    g_domain = {};
    g_lambda = {};
    g_eta = {};
else
    g_domain = domain_generators(x_model, cfg, scales);
    g_lambda = lift_ball_generator(lambda_model, lambda, bounds.R_lambda, scales.lambda);
    g_eta = lift_ball_generator(eta_model, eta, bounds.R_eta, scales.eta);
end

target_pos = V - cfg.epsilon * (x' * x);
[F, params, poly_meta] = add_sos_certificate( ...
    F, params, poly_meta, target_pos, z_current, ...
    [g_cur, g_domain, g_lambda], h_cur, hdeg_cur, dV, cfg, 'current_pos');

target_dec = V - replace(V, x_model, xplus_model) - cfg.alpha * (x' * x);
[F, params, poly_meta] = add_sos_certificate( ...
    F, params, poly_meta, target_dec, z_lift, ...
    [g_cur, g_next, g_domain, g_lambda, g_eta], ...
    [h_cur, h_next], [hdeg_cur, hdeg_next], dV, cfg, 'transition_dec');
end


function [F, params, poly_meta] = add_sos_certificate( ...
    F, params, poly_meta, target_poly, vars, generators, equalities, equality_degrees, dV, cfg, base_name)

cert_sum = 0;
max_term_degree = dV;

for ii = 1:numel(generators)
    [s, c] = polynomial(vars, cfg.sos_multiplier_degree);
    params = [params; c]; %#ok<AGROW>
    sos_constraint = sos(s);
    F = [F, sos_constraint]; %#ok<AGROW>
    poly_meta = append_sos_constraint(poly_meta, sos_constraint, s, vars, ...
        cfg.sos_multiplier_degree, sprintf('%s_g_%02d', base_name, ii), "multiplier");
    cert_sum = cert_sum + s * generators{ii};
    max_term_degree = max(max_term_degree, cfg.sos_multiplier_degree + generator_degree_hint(generators{ii}));
end

for ii = 1:numel(equalities)
    h = equalities{ii};
    hdeg = equality_degrees(ii);
    if hdeg <= 1
        pdeg = stationarity_multiplier_degree(dV, cfg);
    else
        pdeg = cfg.complementarity_multiplier_degree;
    end
    [p, c] = polynomial(vars, pdeg);
    params = [params; c]; %#ok<AGROW>
    F = [F, c <= cfg.coefficient_bound, c >= -cfg.coefficient_bound]; %#ok<AGROW>
    cert_sum = cert_sum + p * h;
    max_term_degree = max(max_term_degree, pdeg + hdeg);
end

residual = target_poly - cert_sum;
residual_degree = next_even(max_term_degree);
sos_constraint = sos(residual);
F = [F, sos_constraint]; %#ok<AGROW>
poly_meta = append_sos_constraint(poly_meta, sos_constraint, residual, vars, ...
    residual_degree, sprintf('%s_residual', base_name), "residual");
end


function pdeg = stationarity_multiplier_degree(dV, cfg)
if dV <= 2
    pdeg = cfg.stationarity_multiplier_degree_dV2;
else
    pdeg = cfg.stationarity_multiplier_degree_dV4;
end
end


function d = generator_degree_hint(poly)
% Conservative hint for this script: all generators are affine or quadratic.
try
    deg = degree(poly);
    if isempty(deg) || isnan(deg)
        d = 2;
    else
        d = max(1, deg);
    end
catch
    d = 2;
end
end


function n = next_even(n)
n = ceil(n);
if mod(n, 2) ~= 0
    n = n + 1;
end
end


function [h, hdeg, g] = kkt_polynomials(xparam, U, lambda_ub, lambda_lb, ...
    lambda_ub_generator, lambda_lb_generator, mpc, cfg, complementarity_scale)
stat = mpc.H * U + mpc.F * xparam + lambda_ub - lambda_lb;
switch string(cfg.stationarity_normalization)
    case "none"
        stat_scale = 1;
    case "H_sigma"
        stat_scale = mpc.H_sigma_max;
    otherwise
        error('Unknown stationarity_normalization: %s', cfg.stationarity_normalization);
end
comp_scale = max(1, complementarity_scale);
h = {
    stat(1) / stat_scale, stat(2) / stat_scale, ...
    lambda_ub(1) * (1 - U(1)) / comp_scale, ...
    lambda_ub(2) * (1 - U(2)) / comp_scale, ...
    lambda_lb(1) * (U(1) + 1) / comp_scale, ...
    lambda_lb(2) * (U(2) + 1) / comp_scale
    };
hdeg = [1, 1, 2, 2, 2, 2];
g = {
    1 - U(1), 1 - U(2), U(1) + 1, U(2) + 1, ...
    lambda_ub_generator(1), lambda_ub_generator(2), ...
    lambda_lb_generator(1), lambda_lb_generator(2)
    };
end


function g = domain_generators(x_model, cfg, scales)
if strcmp(string(scales.mode), "compact_unit_box")
    g = {
        x_model(1) + 1, 1 - x_model(1), ...
        x_model(2) + 1, 1 - x_model(2)
        };
else
    x1_min = cfg.domain_x1(1) / scales.x(1);
    x1_max = cfg.domain_x1(2) / scales.x(1);
    x2_min = cfg.domain_x2(1) / scales.x(2);
    x2_max = cfg.domain_x2(2) / scales.x(2);
    g = {
        x_model(1) - x1_min, x1_max - x_model(1), ...
        x_model(2) - x2_min, x2_max - x_model(2)
        };
end
end


function g = lift_ball_generator(model_var, physical_var, physical_radius, scale_factor)
if scale_factor == 1
    g = {physical_radius^2 - physical_var' * physical_var};
else
    g = {1 - model_var' * model_var};
end
end


function poly_meta = append_sos_constraint(poly_meta, sos_constraint, poly, vars, degree_hint, name, kind)
entry = struct();
entry.name = string(name);
entry.kind = string(kind);
entry.poly = poly;
entry.vars = vars;
entry.degree = degree_hint;
basis = monolist(vars, degree_hint / 2);
entry.gram_dim = length(basis);
entry.is_sos_constraint = true;
entry.sos_constraint = sos_constraint;
entry.gram_source = "sosd_bound_to_constraint_sosid";
entry.candidate_monomial_count = length(basis);
entry.gram_basis_mismatch = false;
poly_meta.sos_meta = append_sos_meta_kkt(poly_meta.sos_meta, entry);
poly_meta.sos_candidate_monomials{end+1} = basis;
end


function sos_meta = empty_sos_meta_kkt()
sos_meta = struct('name', {}, 'kind', {}, 'poly', {}, 'vars', {}, 'degree', {}, ...
    'gram_dim', {}, 'is_sos_constraint', {}, 'sos_constraint', {}, ...
    'gram_source', {}, 'candidate_monomial_count', {}, 'gram_basis_mismatch', {});
end


function meta = append_sos_meta_kkt(meta, entry)
if isempty(meta)
    meta = entry;
else
    meta(end+1) = entry; %#ok<AGROW>
end
end


function sos_oracle = independent_sos_postcheck_kkt(poly_meta, cfg, solver_return_block_count, ...
    solver_monomials, solver_grams, solver_residuals)
if nargin < 3
    solver_return_block_count = NaN;
end
if nargin < 4
    solver_monomials = {};
end
if nargin < 5
    solver_grams = {};
end
if nargin < 6
    solver_residuals = [];
end
sos_oracle = empty_sos_oracle_postcheck();
sos_oracle.method = "meta_bound_gram_min_eig_plus_independent_coeff_match";
sos_oracle.gram_source = "sosd_bound_to_constraint_sosid";
sos_oracle.solver_return_block_count = solver_return_block_count;

sos_meta = poly_meta.sos_meta;
min_gram_eig = inf;
max_abs_coeff_residual = 0;
worst_gram_block = "";
worst_coeff_block = "";
gram_handle_missing_count = 0;
gram_dim_mismatch_count = 0;
block_diagnostics = empty_coeff_debug_row();

for kk = 1:numel(sos_meta)
    [Qv, gram_poly, missing, dim_mismatch] = gram_block_from_meta(sos_meta(kk));
    gram_handle_missing_count = gram_handle_missing_count + double(missing);
    gram_dim_mismatch_count = gram_dim_mismatch_count + double(dim_mismatch);

    block_min_eig = min_gram_eig_from_matrix(Qv);
    if isnan(block_min_eig) || block_min_eig < min_gram_eig
        min_gram_eig = block_min_eig;
        worst_gram_block = sos_meta(kk).name;
    end

    stats = coeff_match_stats(sos_meta(kk).poly, gram_poly, sos_meta(kk).vars);
    block_max_abs_coeff = stats.abs_residual;
    if isnan(block_max_abs_coeff) || block_max_abs_coeff > max_abs_coeff_residual
        max_abs_coeff_residual = block_max_abs_coeff;
        worst_coeff_block = sos_meta(kk).name;
    end

    returned = returned_gram_diagnostics(sos_meta(kk), kk, solver_monomials, solver_grams, solver_residuals);
    block_diag = make_block_diagnostic(sos_meta(kk), kk, Qv, missing, dim_mismatch, ...
        block_min_eig, stats, returned, cfg);
    block_diagnostics = append_debug_struct(block_diagnostics, block_diag);
end

if isinf(min_gram_eig)
    min_gram_eig = NaN;
end

sos_oracle.min_gram_eig = min_gram_eig;
sos_oracle.max_abs_coeff_residual = max_abs_coeff_residual;
sos_oracle.worst_gram_block = worst_gram_block;
sos_oracle.worst_coeff_block = worst_coeff_block;
sos_oracle.gram_handle_missing_count = gram_handle_missing_count;
sos_oracle.gram_dim_mismatch_count = gram_dim_mismatch_count;
sos_oracle.block_diagnostics = block_diagnostics;
sos_oracle.strict_pass = gram_handle_missing_count == 0 && gram_dim_mismatch_count == 0 && ...
    isfinite(min_gram_eig) && isfinite(max_abs_coeff_residual) && ...
    min_gram_eig >= -cfg.gram_eig_tol && max_abs_coeff_residual <= cfg.coeff_match_tol;
sos_oracle.tol_pass = gram_handle_missing_count == 0 && gram_dim_mismatch_count == 0 && ...
    isfinite(min_gram_eig) && isfinite(max_abs_coeff_residual) && ...
    min_gram_eig >= -cfg.sos_oracle_tol_factor * cfg.gram_eig_tol && ...
    max_abs_coeff_residual <= cfg.sos_oracle_tol_factor * cfg.coeff_match_tol;
end


function [Qv, gram_poly, missing, dim_mismatch] = gram_block_from_meta(meta)
missing = isempty(meta.sos_constraint);
dim_mismatch = true;
if missing
    Qv = NaN;
    gram_poly = NaN;
    return;
end
try
    [~, ~, Qraw, Z] = sosd(meta.sos_constraint);
catch
    Qv = NaN;
    gram_poly = NaN;
    missing = true;
    return;
end
Qv = Qraw;
Z = Z(:);
dim_mismatch = isempty(Qv) || size(Qv, 1) ~= length(Z) || size(Qv, 2) ~= length(Z) || ...
    length(Z) ~= meta.gram_dim;
if dim_mismatch || any(isnan(Qv(:)))
    gram_poly = NaN;
    return;
end
gram_poly = Z' * Qv * Z;
end


function block_min_eig = min_gram_eig_from_matrix(Qv)
if isempty(Qv) || any(isnan(Qv(:)))
    block_min_eig = NaN;
elseif isscalar(Qv)
    block_min_eig = Qv;
else
    Qsym = 0.5 * (Qv + Qv');
    block_min_eig = min(eig(Qsym));
end
end


function stats = coeff_match_stats(poly, gram_poly, vars)
stats = struct();
stats.abs_residual = NaN;
stats.rel_residual = NaN;
stats.target_coeff_max_abs = max_abs_coeff(poly, vars);
stats.gram_coeff_max_abs = NaN;
if isnumeric(gram_poly) && any(isnan(gram_poly(:)))
    return;
end
coeff_expr = coefficients(poly - gram_poly, vars);
if isempty(coeff_expr)
    coeff_vals = 0;
else
    coeff_vals = value(coeff_expr);
end
if isempty(coeff_vals)
    coeff_vals = 0;
end
if any(isnan(coeff_vals(:)))
    stats.abs_residual = NaN;
else
    stats.abs_residual = max(abs(coeff_vals(:)));
end
stats.gram_coeff_max_abs = max_abs_coeff(gram_poly, vars);
den = max([1, stats.target_coeff_max_abs, stats.gram_coeff_max_abs]);
stats.rel_residual = stats.abs_residual / den;
end


function max_abs = max_abs_coeff(poly, vars)
if isnumeric(poly) && any(isnan(poly(:)))
    max_abs = NaN;
    return;
end
coeff_expr = coefficients(poly, vars);
if isempty(coeff_expr)
    coeff_vals = 0;
else
    coeff_vals = value(coeff_expr);
end
if isempty(coeff_vals)
    coeff_vals = 0;
end
if any(isnan(coeff_vals(:)))
    max_abs = NaN;
else
    max_abs = max(abs(coeff_vals(:)));
end
end


function returned = returned_gram_diagnostics(meta, block_index, solver_monomials, solver_grams, solver_residuals)
returned = struct();
returned.coeff_match_abs = NaN;
returned.coeff_match_rel = NaN;
returned.gram_min_eig = NaN;
returned.dim_mismatch = true;
returned.reported_residual = NaN;
if block_index > numel(solver_monomials) || block_index > numel(solver_grams)
    return;
end
try
    Z = solver_monomials{block_index};
    Q = solver_grams{block_index};
    if isempty(Z) || isempty(Q)
        return;
    end
    Z = Z(:);
    returned.dim_mismatch = size(Q, 1) ~= length(Z) || size(Q, 2) ~= length(Z) || ...
        length(Z) ~= meta.gram_dim;
    returned.gram_min_eig = min_gram_eig_from_matrix(Q);
    if ~returned.dim_mismatch
        gram_poly = Z' * Q * Z;
        stats = coeff_match_stats(meta.poly, gram_poly, meta.vars);
        returned.coeff_match_abs = stats.abs_residual;
        returned.coeff_match_rel = stats.rel_residual;
    end
catch
    returned.dim_mismatch = true;
end
if isnumeric(solver_residuals) && block_index <= numel(solver_residuals)
    returned.reported_residual = solver_residuals(block_index);
elseif iscell(solver_residuals) && block_index <= numel(solver_residuals)
    try
        returned.reported_residual = solver_residuals{block_index};
    catch
        returned.reported_residual = NaN;
    end
end
end


function block_diag = make_block_diagnostic(meta, block_index, Qv, missing, dim_mismatch, ...
    gram_min_eig, stats, returned, cfg)
block_diag = struct();
block_diag.block_index = block_index;
block_diag.block_name = meta.name;
block_diag.block_kind = meta.kind;
block_diag.block_degree = meta.degree;
block_diag.gram_dim = meta.gram_dim;
block_diag.candidate_monomial_count = meta.candidate_monomial_count;
block_diag.sosd_coeff_match_abs = stats.abs_residual;
block_diag.sosd_coeff_match_rel = stats.rel_residual;
block_diag.target_coeff_max_abs = stats.target_coeff_max_abs;
block_diag.gram_coeff_max_abs = stats.gram_coeff_max_abs;
block_diag.sosd_gram_min_eig = gram_min_eig;
block_diag.sosd_missing = missing;
block_diag.sosd_dim_mismatch = dim_mismatch;
block_diag.returned_coeff_match_abs = returned.coeff_match_abs;
block_diag.returned_coeff_match_rel = returned.coeff_match_rel;
block_diag.returned_gram_min_eig = returned.gram_min_eig;
block_diag.returned_dim_mismatch = returned.dim_mismatch;
block_diag.solvesos_reported_residual = returned.reported_residual;
block_diag.yalmip_check_value = check_single_constraint(meta.sos_constraint);
block_diag.random_eval_max_abs = random_residual_eval(meta.poly, Qv, meta, cfg);
end


function val = check_single_constraint(con)
try
    check_vals = check(con);
    if isempty(check_vals)
        val = NaN;
    else
        val = min(check_vals);
    end
catch
    val = NaN;
end
end


function max_abs_val = random_residual_eval(poly, Qv, meta, cfg)
max_abs_val = NaN;
if isempty(Qv) || any(isnan(Qv(:))) || isempty(meta.sos_constraint)
    return;
end
try
    [~, ~, ~, Z] = sosd(meta.sos_constraint);
    Z = Z(:);
    if size(Qv, 1) ~= length(Z) || size(Qv, 2) ~= length(Z)
        return;
    end
    residual = poly - Z' * Qv * Z;
    rng(11 + meta.gram_dim + meta.degree);
    max_abs_val = 0;
    for ii = 1:cfg.coeff_debug_random_points
        point = -1 + 2 * rand(length(meta.vars), 1);
        val = value(replace(residual, meta.vars, point));
        max_abs_val = max(max_abs_val, abs(val));
    end
catch
    max_abs_val = NaN;
end
end


function rows = empty_coeff_debug_row()
rows = struct('block_index', {}, 'block_name', {}, 'block_kind', {}, ...
    'block_degree', {}, 'gram_dim', {}, 'candidate_monomial_count', {}, ...
    'sosd_coeff_match_abs', {}, 'sosd_coeff_match_rel', {}, ...
    'target_coeff_max_abs', {}, 'gram_coeff_max_abs', {}, ...
    'sosd_gram_min_eig', {}, 'sosd_missing', {}, 'sosd_dim_mismatch', {}, ...
    'returned_coeff_match_abs', {}, 'returned_coeff_match_rel', {}, ...
    'returned_gram_min_eig', {}, 'returned_dim_mismatch', {}, ...
    'solvesos_reported_residual', {}, 'yalmip_check_value', {}, ...
    'random_eval_max_abs', {});
end


function rows = append_debug_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end


function V_fixed = fixed_polynomial(coeffs, monomials)
cv = value(coeffs);
V_fixed = 0;
for ii = 1:numel(cv)
    V_fixed = V_fixed + cv(ii) * monomials(ii);
end
end


function collision = qp_collision_postcheck(V_fixed, xvar, plant, mpc, cfg, scales)
points = collision_points(cfg);
min_decay = inf;
failed_solves = 0;
xp_out_box_count = 0;
worst_x = [NaN; NaN];
worst_u1 = NaN;
worst_xp = [NaN; NaN];
worst_Vx = NaN;
worst_Vxp = NaN;
worst_xp_in_box = false;
for kk = 1:size(points, 2)
    x = points(:, kk);
    [U, ok] = solve_box_qp(x, mpc);
    if ~ok
        failed_solves = failed_solves + 1;
        continue;
    end
    xp = plant.As * x + plant.Bs * U(1);
    Vx = value(replace(V_fixed, xvar, physical_x_to_model_x(x, scales)));
    Vxp = value(replace(V_fixed, xvar, physical_x_to_model_x(xp, scales)));
    decay = Vx - Vxp - cfg.alpha * (x' * x);
    xp_in_box = point_in_domain(xp, cfg);
    xp_out_box_count = xp_out_box_count + double(~xp_in_box);
    if decay < min_decay
        min_decay = decay;
        worst_x = x;
        worst_u1 = U(1);
        worst_xp = xp;
        worst_Vx = Vx;
        worst_Vxp = Vxp;
        worst_xp_in_box = xp_in_box;
    end
end
if isinf(min_decay)
    min_decay = NaN;
end
collision = struct();
collision.point_count = size(points, 2);
collision.failed_solves = failed_solves;
collision.min_decay = min_decay;
collision.xp_out_box_count = xp_out_box_count;
collision.worst_x = worst_x;
collision.worst_u1 = worst_u1;
collision.worst_xp = worst_xp;
collision.worst_Vx = worst_Vx;
collision.worst_Vxp = worst_Vxp;
collision.worst_xp_in_box = worst_xp_in_box;
collision.strict_pass = failed_solves == 0 && isfinite(min_decay) && min_decay >= -cfg.qp_collision_tol;
collision.tol_pass = failed_solves == 0 && isfinite(min_decay) && ...
    min_decay >= -cfg.qp_collision_tol_factor * cfg.qp_collision_tol;
end


function points = collision_points(cfg)
[X1, X2] = meshgrid(cfg.collision_grid_x1, cfg.collision_grid_x2);
points = [X1(:)'; X2(:)'];
corners = [
    cfg.domain_x1(1), cfg.domain_x1(1), cfg.domain_x1(2), cfg.domain_x1(2);
    cfg.domain_x2(1), cfg.domain_x2(2), cfg.domain_x2(1), cfg.domain_x2(2)
    ];
near_origin = cfg.origin_probe_radius * [0 1 -1 0 0; 0 0 0 1 -1];
rng(cfg.collision_random_seed);
rand_points = [
    cfg.domain_x1(1) + diff(cfg.domain_x1) * rand(1, cfg.collision_random_points);
    cfg.domain_x2(1) + diff(cfg.domain_x2) * rand(1, cfg.collision_random_points)
    ];
points = [points, corners, near_origin, rand_points];
end


function collision = qp_collision_roa_postcheck(V_fixed, xvar, plant, mpc, cfg, scales, gamma)
if ~isfinite(gamma)
    collision = struct('point_count', 0, 'failed_solves', 0, 'min_decay', NaN, ...
        'strict_pass', false, 'tol_pass', false);
    return;
end
points = collision_points(cfg);
point_count = 0;
failed_solves = 0;
min_decay = inf;
for kk = 1:size(points, 2)
    x = points(:, kk);
    Vx = value(replace(V_fixed, xvar, physical_x_to_model_x(x, scales)));
    if Vx > gamma + 1e-9
        continue;
    end
    point_count = point_count + 1;
    [U, ok] = solve_box_qp(x, mpc);
    if ~ok
        failed_solves = failed_solves + 1;
        continue;
    end
    xp = plant.As * x + plant.Bs * U(1);
    Vxp = value(replace(V_fixed, xvar, physical_x_to_model_x(xp, scales)));
    decay = Vx - Vxp - cfg.alpha * (x' * x);
    min_decay = min(min_decay, decay);
end
if isinf(min_decay)
    min_decay = NaN;
end
collision = struct();
collision.point_count = point_count;
collision.failed_solves = failed_solves;
collision.min_decay = min_decay;
collision.strict_pass = point_count > 0 && failed_solves == 0 && ...
    isfinite(min_decay) && min_decay >= -cfg.qp_collision_tol;
collision.tol_pass = point_count > 0 && failed_solves == 0 && ...
    isfinite(min_decay) && min_decay >= -cfg.qp_collision_tol_factor * cfg.qp_collision_tol;
end


function x_model = physical_x_to_model_x(x_physical, scales)
x_model = x_physical(:) ./ scales.x(:);
end


function inside = point_in_domain(x, cfg)
inside = x(1) >= cfg.domain_x1(1) - 1e-12 && x(1) <= cfg.domain_x1(2) + 1e-12 && ...
    x(2) >= cfg.domain_x2(1) - 1e-12 && x(2) <= cfg.domain_x2(2) + 1e-12;
end


function [U, ok] = solve_box_qp(x, mpc)
ok = true;
try
    opts = optimoptions('quadprog', 'Display', 'off');
    [U, ~, exitflag] = quadprog(mpc.H, mpc.F * x, [], [], [], [], ...
        mpc.u_min * ones(mpc.N, 1), mpc.u_max * ones(mpc.N, 1), [], opts);
catch
    try
        [U, ~, exitflag] = quadprog(mpc.H, mpc.F * x, [], [], [], [], ...
            mpc.u_min * ones(mpc.N, 1), mpc.u_max * ones(mpc.N, 1));
    catch
        U = NaN(mpc.N, 1);
        ok = false;
        return;
    end
end
ok = ok && ~isempty(U) && exitflag > 0;
end


function roa = find_roa_gamma(V_fixed, xvar, cfg, scales)
roa = struct('gamma', NaN, 'status', "not_attempted", ...
    'strict_pass', false, 'tol_pass', false);
boundary_points = boundary_sample_points(cfg);
Vb = zeros(1, size(boundary_points, 2));
for kk = 1:size(boundary_points, 2)
    Vb(kk) = value(replace(V_fixed, xvar, physical_x_to_model_x(boundary_points(:, kk), scales)));
end
hi = 0.95 * min(Vb);
if ~isfinite(hi) || hi <= cfg.roa_gamma_floor
    roa.status = "no_positive_boundary_level";
    return;
end
lo = 0;
best = 0;
for kk = 1:cfg.roa_bisection_steps
    mid = 0.5 * (lo + hi);
    if roa_inclusion_sos_feasible(V_fixed, xvar, mid, cfg, scales)
        best = mid;
        lo = mid;
    else
        hi = mid;
    end
end
if best > cfg.roa_gamma_floor
    roa.gamma = best;
    roa.status = "verified_roa_subset_by_sos";
    roa.strict_pass = true;
    roa.tol_pass = true;
else
    roa.status = "roa_inclusion_not_verified";
end
end


function points = boundary_sample_points(cfg)
x1 = linspace(cfg.domain_x1(1), cfg.domain_x1(2), 41);
x2 = linspace(cfg.domain_x2(1), cfg.domain_x2(2), 41);
top = [x1; cfg.domain_x2(2) * ones(size(x1))];
bottom = [x1; cfg.domain_x2(1) * ones(size(x1))];
left = [cfg.domain_x1(1) * ones(size(x2)); x2];
right = [cfg.domain_x1(2) * ones(size(x2)); x2];
points = [top, bottom, left, right];
end


function feasible = roa_inclusion_sos_feasible(V_fixed, xvar, gamma, cfg, scales)
F = [];
params = [];
candidates = {};
psi = domain_generators(xvar, cfg, scales);
for ii = 1:numel(psi)
    [s, c] = polynomial(xvar, cfg.domain_multiplier_degree);
    params = [params; c]; %#ok<AGROW>
    scon = sos(s);
    rcon = sos(psi{ii} - s * (gamma - V_fixed));
    F = [F, scon, rcon]; %#ok<AGROW>
    candidates{end+1} = monolist(xvar, cfg.domain_multiplier_degree / 2); %#ok<AGROW>
    candidates{end+1} = monolist(xvar, next_even(cfg.domain_multiplier_degree + max(2, cfg.domain_multiplier_degree)) / 2); %#ok<AGROW>
end
ops = make_solver_options(cfg);
ops.verbose = 0;
try
    sol = solvesos(F, [], ops, params, candidates);
    feasible = sol.problem == 0;
catch
    feasible = false;
end
end


function plant = sampled_execution_model(Ts)
plant = struct();
plant.Ts = Ts;
plant.As = [1, Ts; 0, 1];
plant.Bs = [0.5 * Ts^2; Ts];
end


function bounds = build_lift_bounds(plant, mpc, cfg)
if cfg.global_attempt_mode
    bounds = struct();
    bounds.Rx = Inf;
    bounds.Ru = sqrt(mpc.N);
    bounds.R_lambda = Inf;
    bounds.R_eta = Inf;
    return;
end
Rx = norm([max(abs(cfg.domain_x1)); max(abs(cfg.domain_x2))], 2);
Ru = sqrt(mpc.N);
Rlambda = norm(mpc.H, 2) * Ru + norm(mpc.F, 2) * Rx;
Rx_next = norm(plant.As, 2) * Rx + norm(plant.Bs, 2) * 1;
Reta = norm(mpc.H, 2) * Ru + norm(mpc.F, 2) * Rx_next;
bounds = struct();
bounds.Rx = Rx;
bounds.Ru = Ru;
bounds.R_lambda = 1.05 * Rlambda;
bounds.R_eta = 1.05 * Reta;
end


function scales = build_variable_scaling(bounds, cfg)
scales = struct();
scales.mode = string(cfg.lift_variable_scaling);
switch scales.mode
    case "none"
        scales.x = [1; 1];
        scales.U = ones(2, 1);
        scales.lambda = 1;
        scales.eta = 1;
    case "compact_unit_box"
        scales.x = [max(abs(cfg.domain_x1)); max(abs(cfg.domain_x2))];
        scales.U = ones(2, 1);
        scales.lambda = bounds.R_lambda;
        scales.eta = bounds.R_eta;
    case "fixed_state_20_5"
        scales.x = cfg.fixed_state_scale(:);
        scales.U = ones(2, 1);
        scales.lambda = 1;
        scales.eta = 1;
    otherwise
        error('Unknown lift_variable_scaling: %s', cfg.lift_variable_scaling);
end
end


function mpc = build_standard_mpc_data()
mpc.Ad = [1 1; 0 1];
mpc.Bd = [0; 1];
mpc.C = [1 0];
mpc.Q = mpc.C' * mpc.C;
mpc.R = 0.1;
mpc.N = 2;
mpc.u_min = -1;
mpc.u_max = 1;
if exist('dare', 'file') == 2
    mpc.P = dare(mpc.Ad, mpc.Bd, mpc.Q, mpc.R);
else
    mpc.P = solve_dare_fixed_point(mpc.Ad, mpc.Bd, mpc.Q, mpc.R);
end
mpc.P = 0.5 * (mpc.P + mpc.P');
[mpc.H, mpc.F] = build_condensed_qp_matrices(mpc);
mpc.H = 0.5 * (mpc.H + mpc.H');
mpc.H_sigma_max = max(svd(mpc.H));
mpc.H_normalized = mpc.H / mpc.H_sigma_max;
mpc.F_normalized = mpc.F / mpc.H_sigma_max;
end


function P = solve_dare_fixed_point(A, B, Q, R)
P = Q;
for iter = 1:10000
    P_next = A' * P * A - A' * P * B * ((R + B' * P * B) \ (B' * P * A)) + Q;
    P_next = 0.5 * (P_next + P_next');
    if norm(P_next - P, 'fro') <= 1e-13
        P = P_next;
        return;
    end
    P = P_next;
end
error('Fixed-point DARE iteration did not converge.');
end


function [H, F] = build_condensed_qp_matrices(mpc)
n = size(mpc.Ad, 1);
N = mpc.N;
H = zeros(N, N);
F = zeros(N, n);
for k = 0:N-1
    [Phi_k, Gamma_k] = prediction_matrices(mpc.Ad, mpc.Bd, N, k);
    H = H + 2 * (Gamma_k' * mpc.Q * Gamma_k);
    F = F + 2 * (Gamma_k' * mpc.Q * Phi_k);
    H(k+1,k+1) = H(k+1,k+1) + 2 * mpc.R;
end
[Phi_N, Gamma_N] = prediction_matrices(mpc.Ad, mpc.Bd, N, N);
H = H + 2 * (Gamma_N' * mpc.P * Gamma_N);
F = F + 2 * (Gamma_N' * mpc.P * Phi_N);
end


function [Phi, Gamma] = prediction_matrices(A, B, N, k)
n = size(A, 1);
Phi = A^k;
Gamma = zeros(n, N);
for j = 1:min(k, N)
    Gamma(:,j) = A^(k-j) * B;
end
end


function precheck = precheck_standard_mpc_data(mpc, cfg)
precheck = struct();
precheck.H_min_eig = min(eig(0.5 * (mpc.H + mpc.H')));
precheck.box_nonempty = all(mpc.u_min < mpc.u_max);
precheck.origin_qp_equilibrium = norm(mpc.F * [0; 0]) <= 1e-12;
tool_names = {'sdpvar','sos','solvesos','polynomial','coefficients','quadprog'};
missing = strings(0, 1);
for kk = 1:numel(tool_names)
    if isempty(which(tool_names{kk}))
        missing(end+1, 1) = string(tool_names{kk}); %#ok<AGROW>
    end
end
precheck.tools_available = isempty(missing);
precheck.tool_info = strjoin(missing, ",");
precheck.cfg_solver = cfg.solver;
end


function row = base_row(Ts, dV, cfg, precheck)
row = struct();
row.Ts = Ts;
row.d_V = dV;
row.epsilon = cfg.epsilon;
row.alpha = cfg.alpha;
row.domain_x1_min = cfg.domain_x1(1);
row.domain_x1_max = cfg.domain_x1(2);
row.domain_x2_min = cfg.domain_x2(1);
row.domain_x2_max = cfg.domain_x2(2);
row.kkt_mode = cfg.kkt_mode;
row.relaxed_delta = cfg.relaxed_delta;
row.certificate_claim_mode = cfg.certificate_claim_mode;
row.require_invariance_for_status = cfg.require_invariance_for_status;
row.global_attempt_mode = cfg.global_attempt_mode;
row.noncompact_lifted_set = cfg.noncompact_lifted_set;
row.analysis_domain = cfg.analysis_domain;
row.solver = cfg.solver;
row.run_label = cfg.run_label;
row.domain_profile = cfg.domain_profile;
row.global_sanity_x1_min = cfg.global_sanity_x1(1);
row.global_sanity_x1_max = cfg.global_sanity_x1(2);
row.global_sanity_x2_min = cfg.global_sanity_x2(1);
row.global_sanity_x2_max = cfg.global_sanity_x2(2);
row.yalmip_sos_model = cfg.yalmip_sos_model;
row.yalmip_sos_scale = cfg.yalmip_sos_scale;
row.stationarity_normalization = cfg.stationarity_normalization;
row.lift_variable_scaling = cfg.lift_variable_scaling;
row.fixed_state_scale_x1 = cfg.fixed_state_scale(1);
row.fixed_state_scale_x2 = cfg.fixed_state_scale(2);
row.H_min_eig = precheck.H_min_eig;
row.box_nonempty = precheck.box_nonempty;
row.origin_qp_equilibrium = precheck.origin_qp_equilibrium;
row.certificate_status = "not_run";
row.failure_stage = "";
row.solver_problem_code = NaN;
row.solver_info = "";
row.execution_prediction_mismatch_norm = NaN;
row.R_lambda = NaN;
row.R_eta = NaN;
row.lifted_variable_count = NaN;
row.sos_constraint_count = 0;
row.total_gram_dim_estimate = NaN;
row.yalmip_check_min = NaN;
row.yalmip_check_strict_pass = false;
row.yalmip_check_tol_pass = false;
row.gram_min_eig = NaN;
row.gram_worst_block = "";
row.coeff_match_max_abs = NaN;
row.coeff_match_worst_block = "";
row.sos_postcheck_gram_handle_missing_count = NaN;
row.sos_postcheck_gram_dim_mismatch_count = NaN;
row.sos_oracle_strict_pass = false;
row.sos_oracle_tol_pass = false;
row.qp_collision_point_count = 0;
row.qp_collision_min_decay = NaN;
row.qp_collision_strict_pass = false;
row.qp_collision_tol_pass = false;
row.qp_collision_xp_out_box_count = NaN;
row.qp_collision_worst_x1 = NaN;
row.qp_collision_worst_x2 = NaN;
row.qp_collision_worst_u1 = NaN;
row.qp_collision_worst_xp1 = NaN;
row.qp_collision_worst_xp2 = NaN;
row.qp_collision_worst_Vx = NaN;
row.qp_collision_worst_Vxp = NaN;
row.qp_collision_worst_xp_in_box = false;
row.qp_collision_roa_point_count = 0;
row.qp_collision_roa_min_decay = NaN;
row.qp_collision_roa_strict_pass = false;
row.qp_collision_roa_tol_pass = false;
row.roa_gamma = NaN;
row.roa_invariance_status = "";
row.roa_inclusion_strict_pass = false;
row.roa_inclusion_tol_pass = false;
row.solve_time_sec = NaN;
end


function cert = empty_cert(Ts, dV, cfg)
cert = struct();
cert.Ts = Ts;
cert.d_V = dV;
cert.cfg_run_label = cfg.run_label;
cert.status = "not_run";
end


function rows = append_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end


function rows = append_coeff_debug_rows(rows, cert, result_row, cfg)
if ~cfg.coeff_debug_enabled || ~isfield(cert, 'sos_oracle') || ...
        ~isfield(cert.sos_oracle, 'block_diagnostics')
    return;
end
diag = cert.sos_oracle.block_diagnostics;
for kk = 1:numel(diag)
    debug_row = struct();
    debug_row.Ts = result_row.Ts;
    debug_row.d_V = result_row.d_V;
    debug_row.run_label = result_row.run_label;
    debug_row.domain_profile = result_row.domain_profile;
    debug_row.yalmip_sos_model = cfg.yalmip_sos_model;
    debug_row.yalmip_sos_scale = cfg.yalmip_sos_scale;
    debug_row.stationarity_normalization = cfg.stationarity_normalization;
    debug_row.lift_variable_scaling = cfg.lift_variable_scaling;
    fields = fieldnames(diag(kk));
    for ff = 1:numel(fields)
        debug_row.(fields{ff}) = diag(kk).(fields{ff});
    end
    if isempty(rows)
        rows = debug_row;
    else
        rows(end+1) = debug_row; %#ok<AGROW>
    end
end
end


function print_case_summary(row)
fprintf('  Ts=%.12g dV=%d status=%s, solver=%g, gram_min=%.4g, coeff=%.4g, qp_min=%.4g, gamma=%.4g, time=%.2fs\n', ...
    row.Ts, row.d_V, row.certificate_status, row.solver_problem_code, ...
    row.gram_min_eig, row.coeff_match_max_abs, row.qp_collision_min_decay, ...
    row.roa_gamma, row.solve_time_sec);
end


function meta = strip_runtime_handles(poly_meta)
meta = poly_meta;
if isfield(meta, 'sos_meta')
    for kk = 1:numel(meta.sos_meta)
        meta.sos_meta(kk).sos_constraint = [];
        meta.sos_meta(kk).poly = [];
        meta.sos_meta(kk).vars = [];
    end
end
if isfield(meta, 'sos_candidate_monomials')
    meta.sos_candidate_monomials = {};
end
fields_to_clear = {'x','x_physical','U','lambda','Un','eta','z_current','z_lift','V','V_coeff','V_monomials'};
for kk = 1:numel(fields_to_clear)
    if isfield(meta, fields_to_clear{kk})
        meta.(fields_to_clear{kk}) = [];
    end
end
end


function sos_oracle = empty_sos_oracle_postcheck()
sos_oracle = struct();
sos_oracle.method = "not_run";
sos_oracle.gram_source = "";
sos_oracle.solver_return_block_count = NaN;
sos_oracle.min_gram_eig = NaN;
sos_oracle.max_abs_coeff_residual = NaN;
sos_oracle.worst_gram_block = "";
sos_oracle.worst_coeff_block = "";
sos_oracle.gram_handle_missing_count = 0;
sos_oracle.gram_dim_mismatch_count = 0;
sos_oracle.block_diagnostics = empty_coeff_debug_row();
sos_oracle.strict_pass = false;
sos_oracle.tol_pass = false;
end


function script_dir = get_script_dir()
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
if isempty(script_dir)
    script_dir = pwd;
end
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


function ensure_yalmip_export_precedence()
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
