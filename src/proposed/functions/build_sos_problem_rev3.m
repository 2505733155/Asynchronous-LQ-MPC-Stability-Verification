function [F, params, poly_meta] = build_sos_problem_rev3(reg_poly, Phi, q, transitions, origin_regions, dV, dmult, degree_meta, cfg)
%BUILD_SOS_PROBLEM_REV3 Build the Rev3 piecewise-polynomial SOS model.
% Each SOS block stores the corresponding YALMIP SOS constraint object, so
% postchecks can retrieve its Gram decomposition with sosd(constraint).

x = sdpvar(2, 1);
xnorm2 = x' * x;
M = numel(reg_poly);
F = [];
params = [];
poly_meta = struct();
poly_meta.x = x;
poly_meta.p = cell(1, M);
poly_meta.p_coeff = cell(1, M);
poly_meta.dV = dV;
poly_meta.dmult = dmult;
poly_meta.epsilon = cfg.epsilon;
poly_meta.alpha = cfg.alpha;
poly_meta.sos_meta = empty_sos_meta_rev3();
poly_meta.sos_candidate_monomials = {};

for i = 1:M
    [p, c] = polynomial(x, dV);
    poly_meta.p{i} = p;
    poly_meta.p_coeff{i} = c;
    params = [params; c]; %#ok<AGROW>
    F = [F, c <= cfg.coefficient_bound, c >= -cfg.coefficient_bound]; %#ok<AGROW>
    if ismember(i, origin_regions)
        gradp = jacobian(p, x);
        F = [F, replace(p, x, [0; 0]) == 0, ...
            replace(gradp(1), x, [0; 0]) == 0, ...
            replace(gradp(2), x, [0; 0]) == 0]; %#ok<AGROW>
    end
end

for i = 1:M
    g = poly_generators_rev3(reg_poly(i).H, reg_poly(i).b, x);
    fpos = poly_meta.p{i} - cfg.epsilon * xnorm2;
    [F, params, poly_meta.sos_meta, poly_meta.sos_candidate_monomials] = add_sos_on_generators_rev3( ...
        F, params, poly_meta.sos_meta, poly_meta.sos_candidate_monomials, ...
        fpos, g, x, dmult, degree_meta.d_sigma, sprintf('region_pos_%02d', i));
end

for kk = 1:numel(transitions)
    i = transitions(kk).i;
    j = transitions(kk).j;
    xp = Phi{i} * x + q{i};
    pnext = replace(poly_meta.p{j}, x, xp);
    fdec = poly_meta.p{i} - pnext - cfg.alpha * xnorm2;
    g = poly_generators_rev3(transitions(kk).G, transitions(kk).w, x);
    [F, params, poly_meta.sos_meta, poly_meta.sos_candidate_monomials] = add_sos_on_generators_rev3( ...
        F, params, poly_meta.sos_meta, poly_meta.sos_candidate_monomials, ...
        fdec, g, x, dmult, degree_meta.d_sigma, sprintf('transition_dec_%02d_%02d_%02d', kk, i, j));
end

switch string(cfg.normalization_mode)
    case "fixed_margin_plus_box"
        % Scale is anchored by fixed epsilon/alpha and coefficient box.
    otherwise
        error('Unknown normalization_mode: %s', cfg.normalization_mode);
end

poly_meta.degree_meta = degree_meta;
assert(numel(poly_meta.sos_candidate_monomials) == nnz([poly_meta.sos_meta.is_sos_constraint]), ...
    'Candidate monomial list must match the number of SOS constraints.');
end


function [F, params, sos_meta, candidate_monomials] = add_sos_on_generators_rev3( ...
    F, params, sos_meta, candidate_monomials, target_poly, generators, x, dmult, residual_degree, base_name)

mult_sum = 0;
for ell = 1:numel(generators)
    if dmult == 0
        s = sdpvar(1, 1);
        F = [F, s >= 0]; %#ok<AGROW>
        params = [params; s]; %#ok<AGROW>
        sos_meta = append_scalar_nonnegative_meta_rev3(sos_meta, ...
            sprintf('%s_mult_%02d', base_name, ell), "multiplier_scalar", s);
    else
        [s, c] = polynomial(x, dmult);
        params = [params; c]; %#ok<AGROW>
        basis = monolist(x, dmult / 2);
        sos_constraint = sos(s);
        F = [F, sos_constraint]; %#ok<AGROW>
        meta_entry = make_sos_meta_entry_rev3(sprintf('%s_mult_%02d', base_name, ell), ...
            "multiplier", s, dmult, true);
        meta_entry.sos_constraint = sos_constraint;
        meta_entry.gram_source = "sosd_bound_to_constraint_sosid";
        meta_entry.candidate_monomial_count = length(basis);
        meta_entry.gram_basis_mismatch = meta_entry.candidate_monomial_count ~= meta_entry.gram_dim;
        sos_meta = append_sos_meta_rev3(sos_meta, meta_entry);
        candidate_monomials{end+1} = basis; %#ok<AGROW>
    end
    mult_sum = mult_sum + s * generators{ell};
end

residual = target_poly - mult_sum;
basis = monolist(x, residual_degree / 2);
sos_constraint = sos(residual);
F = [F, sos_constraint]; %#ok<AGROW>
meta_entry = make_sos_meta_entry_rev3(sprintf('%s_residual', base_name), ...
    "residual", residual, residual_degree, true);
meta_entry.sos_constraint = sos_constraint;
meta_entry.gram_source = "sosd_bound_to_constraint_sosid";
meta_entry.candidate_monomial_count = length(basis);
meta_entry.gram_basis_mismatch = meta_entry.candidate_monomial_count ~= meta_entry.gram_dim;
sos_meta = append_sos_meta_rev3(sos_meta, meta_entry);
candidate_monomials{end+1} = basis; %#ok<AGROW>
end


function meta_entry = make_sos_meta_entry_rev3(name, kind, poly, degree, is_sos_constraint)
meta_entry = struct();
meta_entry.name = string(name);
meta_entry.kind = string(kind);
meta_entry.poly = poly;
meta_entry.degree = degree;
meta_entry.gram_dim = monomial_count_2d_rev3(ceil(degree / 2));
meta_entry.is_sos_constraint = is_sos_constraint;
meta_entry.sos_constraint = [];
meta_entry.gram_source = "none";
meta_entry.candidate_monomial_count = 0;
meta_entry.gram_basis_mismatch = true;
end


function meta = empty_sos_meta_rev3()
meta = struct('name', {}, 'kind', {}, 'poly', {}, 'degree', {}, 'gram_dim', {}, ...
    'is_sos_constraint', {}, 'sos_constraint', {}, ...
    'gram_source', {}, 'candidate_monomial_count', {}, 'gram_basis_mismatch', {});
end


function meta = append_sos_meta_rev3(meta, entry)
if isempty(meta)
    meta = entry;
else
    meta(end+1) = entry; %#ok<AGROW>
end
end


function meta = append_scalar_nonnegative_meta_rev3(meta, name, kind, scalar_var)
entry = struct();
entry.name = string(name);
entry.kind = string(kind);
entry.poly = scalar_var;
entry.degree = 0;
entry.gram_dim = 1;
entry.is_sos_constraint = false;
entry.sos_constraint = [];
entry.gram_source = "scalar_nonnegative";
entry.candidate_monomial_count = 1;
entry.gram_basis_mismatch = false;
meta = append_sos_meta_rev3(meta, entry);
end


function g = poly_generators_rev3(G, w, x)
m = size(G, 1);
g = cell(1, m);
for ell = 1:m
    g{ell} = w(ell) - G(ell, :) * x;
end
end


function n = monomial_count_2d_rev3(d)
n = nchoosek(2 + d, d);
end
