function proof = verify_forward_invariance_union_rev3(reg, Phi, q, scale_ctx, cfg)
%VERIFY_FORWARD_INVARIANCE_UNION_REV3 Check PWA domain forward invariance.
% The domain is the region union D = union_i {x | H_i x <= b_i}. If the
% partition is already verified as full-space, forward invariance is
% immediate. Otherwise, this function maximizes the common escape margin
% for each complement selection. Floating-point LP solutions are treated
% as tolerance-qualified evidence unless a validated upper bound is used.

proof = empty_forward_invariance_proof_rev3();
proof.checked = true;

if scale_ctx.partition_full_space_covered
    proof.verified = true;
    proof.failed = false;
    proof.status = "verified_full_state_space_partition";
    proof.method = "full_state_space_partition";
    proof.lp_count = 0;
    return;
end

proof.method = "region_union_max_escape_margin_lp";
M = numel(reg);
n = size(reg(1).H, 2);
tau = cfg.forward_invariance_tol;
opts = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
if isfield(cfg, 'forward_invariance_linprog_fcn') && ~isempty(cfg.forward_invariance_linprog_fcn)
    linprog_fcn = cfg.forward_invariance_linprog_fcn;
else
    linprog_fcn = @linprog;
end

for i = 1:M
    choices = ones(1, M);
    while true
        G = [reg(i).H, zeros(size(reg(i).H, 1), 1)];
        w = reg(i).b;
        for j = 1:M
            ell = choices(j);
            h = reg(j).H(ell, :);
            G = [G; -h * Phi{i}, 1]; %#ok<AGROW>
            w = [w; -reg(j).b(ell) + h * q{i}]; %#ok<AGROW>
        end

        objective = [zeros(n, 1); -1];
        [z_opt, fval, exitflag] = linprog_fcn(objective, G, w, [], [], [], [], opts);
        proof.lp_count = proof.lp_count + 1;
        has_finite_candidate = numel(z_opt) == n + 1 && all(isfinite(z_opt));
        if has_finite_candidate
            x_candidate = z_opt(1:n);
            [valid_counterexample, actual_margin] = ...
                is_valid_forward_invariance_counterexample_rev3( ...
                reg, Phi, q, i, x_candidate, tau);
            if valid_counterexample
                proof.verified = false;
                proof.tolerance_verified = false;
                proof.failed = true;
                proof.status = "failed_forward_invariant_region_union";
                proof.counterexample_region = i;
                proof.counterexample_x = x_candidate;
                proof.counterexample_xp = Phi{i} * x_candidate + q{i};
                proof.counterexample_choice = join(string(choices), ";");
                proof.counterexample_exitflag = exitflag;
                proof.counterexample_escape_margin = actual_margin;
                return;
            end
        end

        if exitflag == 1 && has_finite_candidate && isfinite(fval)
            delta_opt = z_opt(end);
            proof = update_max_escape_margin_rev3(proof, delta_opt, i, choices);
            if delta_opt > tau
                proof.invalid_candidate_count = proof.invalid_candidate_count + 1;
                proof = record_ambiguous_lp_rev3(proof, i, choices, exitflag);
            end
        else
            proof = record_ambiguous_lp_rev3(proof, i, choices, exitflag);
        end

        [choices, done] = next_choice_vector_rev3(choices, reg);
        if done
            break;
        end
    end
end

if proof.ambiguous_lp_count > 0
    proof.verified = false;
    proof.failed = false;
    proof.status = "forward_invariance_not_verified";
    return;
end

proof.verified = false;
proof.tolerance_verified = true;
proof.failed = false;
proof.status = "forward_invariance_verified_with_tolerance";
end


function [tf, min_region_escape_margin] = ...
        is_valid_forward_invariance_counterexample_rev3(reg, Phi, q, i, x_cex, tau)
min_region_escape_margin = NaN;
if isempty(x_cex) || any(~isfinite(x_cex))
    tf = false;
    return;
end

start_slack = reg(i).H * x_cex - reg(i).b;
if any(start_slack > tau)
    tf = false;
    return;
end

xp = Phi{i} * x_cex + q{i};
min_region_escape_margin = inf;
for j = 1:numel(reg)
    region_escape_margin = max(reg(j).H * xp - reg(j).b);
    min_region_escape_margin = min(min_region_escape_margin, region_escape_margin);
end
tf = min_region_escape_margin > tau;
end


function proof = update_max_escape_margin_rev3(proof, delta_opt, i, choices)
if isnan(proof.max_escape_margin) || delta_opt > proof.max_escape_margin
    proof.max_escape_margin = delta_opt;
    proof.max_escape_margin_region = i;
    proof.max_escape_margin_choice = join(string(choices), ";");
end
end


function proof = record_ambiguous_lp_rev3(proof, i, choices, exitflag)
proof.ambiguous_lp_count = proof.ambiguous_lp_count + 1;
if isnan(proof.counterexample_exitflag)
    proof.counterexample_region = i;
    proof.counterexample_choice = join(string(choices), ";");
    proof.counterexample_exitflag = exitflag;
end
end


function proof = empty_forward_invariance_proof_rev3()
proof = struct();
proof.checked = false;
proof.verified = false;
proof.tolerance_verified = false;
proof.failed = false;
proof.status = "not_checked";
proof.method = "not_run";
proof.lp_count = 0;
proof.counterexample_region = NaN;
proof.counterexample_x = [NaN; NaN];
proof.counterexample_xp = [NaN; NaN];
proof.counterexample_choice = "";
proof.counterexample_exitflag = NaN;
proof.counterexample_escape_margin = NaN;
proof.invalid_candidate_count = 0;
proof.ambiguous_lp_count = 0;
proof.max_escape_margin = NaN;
proof.max_escape_margin_region = NaN;
proof.max_escape_margin_choice = "";
end


function [choices, done] = next_choice_vector_rev3(choices, reg)
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
