function sos_oracle = independent_sos_postcheck_rev3(poly_meta, cfg, solver_return_block_count)
%INDEPENDENT_SOS_POSTCHECK_REV3 Check explicit Gram blocks bound in sos_meta.
%
% The check intentionally avoids solvesos-returned gram/monomial cell
% indexing. Each SOS block is reconstructed from the SOS constraint handle
% stored in sos_meta(k), using YALMIP's sosd lookup by sosid.

if nargin < 3
    solver_return_block_count = NaN;
end

sos_oracle = empty_sos_oracle_postcheck_rev3();
sos_oracle.method = "meta_bound_gram_min_eig_plus_independent_coeff_match";
sos_oracle.gram_source = "sosd_bound_to_constraint_sosid";
sos_oracle.solver_return_block_count = solver_return_block_count;

sos_meta = poly_meta.sos_meta;
x = poly_meta.x;
if isempty(sos_meta)
    return;
end

min_gram_eig = inf;
max_abs_coeff_residual = 0;
worst_gram_block = "";
worst_coeff_block = "";
sos_block_count = 0;
gram_handle_missing_count = 0;
gram_dim_mismatch_count = 0;

for kk = 1:numel(sos_meta)
    if sos_meta(kk).is_sos_constraint
        sos_block_count = sos_block_count + 1;
        [Qv, gram_poly, missing, dim_mismatch] = gram_block_from_meta_rev3(sos_meta(kk));
        gram_handle_missing_count = gram_handle_missing_count + double(missing);
        gram_dim_mismatch_count = gram_dim_mismatch_count + double(dim_mismatch);
    else
        Qv = value(sos_meta(kk).poly);
        gram_poly = sos_meta(kk).poly;
    end

    block_min_eig = min_gram_eig_from_matrix_rev3(Qv);
    if isnan(block_min_eig) || block_min_eig < min_gram_eig
        min_gram_eig = block_min_eig;
        worst_gram_block = sos_meta(kk).name;
    end

    if sos_meta(kk).is_sos_constraint
        block_max_abs_coeff = coeff_match_residual_rev3(sos_meta(kk).poly, gram_poly, x);
    else
        block_max_abs_coeff = 0;
    end
    if isnan(block_max_abs_coeff) || block_max_abs_coeff > max_abs_coeff_residual
        max_abs_coeff_residual = block_max_abs_coeff;
        worst_coeff_block = sos_meta(kk).name;
    end
end

if isinf(min_gram_eig)
    min_gram_eig = NaN;
end

sos_oracle.min_gram_eig = min_gram_eig;
sos_oracle.max_abs_coeff_residual = max_abs_coeff_residual;
sos_oracle.worst_gram_block = worst_gram_block;
sos_oracle.worst_coeff_block = worst_coeff_block;
sos_oracle.sos_block_count = sos_block_count;
sos_oracle.gram_handle_missing_count = gram_handle_missing_count;
sos_oracle.gram_dim_mismatch_count = gram_dim_mismatch_count;
sos_oracle.strict_pass = gram_handle_missing_count == 0 && gram_dim_mismatch_count == 0 && ...
    isfinite(min_gram_eig) && isfinite(max_abs_coeff_residual) && ...
    min_gram_eig >= -cfg.gram_eig_tol && max_abs_coeff_residual <= cfg.coeff_match_tol;
sos_oracle.tol_pass = gram_handle_missing_count == 0 && gram_dim_mismatch_count == 0 && ...
    isfinite(min_gram_eig) && isfinite(max_abs_coeff_residual) && ...
    min_gram_eig >= -cfg.sos_oracle_tol_factor * cfg.gram_eig_tol && ...
    max_abs_coeff_residual <= cfg.sos_oracle_tol_factor * cfg.coeff_match_tol;
end


function [Qv, gram_poly, missing, dim_mismatch] = gram_block_from_meta_rev3(meta)
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


function block_min_eig = min_gram_eig_from_matrix_rev3(Qv)
if isempty(Qv) || any(isnan(Qv(:)))
    block_min_eig = NaN;
elseif isscalar(Qv)
    block_min_eig = Qv;
else
    Qsym = 0.5 * (Qv + Qv');
    block_min_eig = min(eig(Qsym));
end
end


function block_max_abs_coeff = coeff_match_residual_rev3(poly, gram_poly, x)
if isnumeric(gram_poly) && any(isnan(gram_poly(:)))
    block_max_abs_coeff = NaN;
    return;
end

coeff_expr = coefficients(poly - gram_poly, x);
if isempty(coeff_expr)
    coeff_vals = 0;
else
    coeff_vals = value(coeff_expr);
end
if isempty(coeff_vals)
    coeff_vals = 0;
end
if any(isnan(coeff_vals(:)))
    block_max_abs_coeff = NaN;
else
    block_max_abs_coeff = max(abs(coeff_vals(:)));
end
end


function sos_oracle = empty_sos_oracle_postcheck_rev3()
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
