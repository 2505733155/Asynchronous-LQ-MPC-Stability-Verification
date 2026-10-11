classdef SosPostcheckGramBindingRev3Test < matlab.unittest.TestCase
    % Regression tests for sosid-bound Gram postchecks.

    methods (Test)
        function solvedMinimalProblemUsesMetaBoundGramHandles(testCase)
            [F, params, polyMeta, cfg] = SosPostcheckGramBindingRev3Test.buildMinimalProblem(4, 2);
            [sol, ~, sosGrams] = solvesos(F, [], SosPostcheckGramBindingRev3Test.options(), ...
                params, polyMeta.sos_candidate_monomials);
            testCase.verifyEqual(sol.problem, 0, sol.info);

            oracle = independent_sos_postcheck_rev3(polyMeta, cfg, numel(sosGrams));
            testCase.verifyEqual(oracle.method, ...
                "meta_bound_gram_min_eig_plus_independent_coeff_match");
            testCase.verifyEqual(oracle.gram_source, "sosd_bound_to_constraint_sosid");
            testCase.verifyEqual(oracle.sos_block_count, ...
                nnz([polyMeta.sos_meta.is_sos_constraint]));
            testCase.verifyEqual(oracle.gram_handle_missing_count, 0);
            testCase.verifyEqual(oracle.gram_dim_mismatch_count, 0);
            testCase.verifyTrue(oracle.tol_pass);
        end

        function sosMetaContainsDirectGramHandles(testCase)
            [~, ~, polyMeta] = SosPostcheckGramBindingRev3Test.buildMinimalProblem(4, 2);
            sosMask = [polyMeta.sos_meta.is_sos_constraint];
            sosMeta = polyMeta.sos_meta(sosMask);

            testCase.verifyTrue(all(arrayfun(@(m) ~isempty(m.sos_constraint), sosMeta)));
            testCase.verifyEqual([sosMeta.candidate_monomial_count], [sosMeta.gram_dim]);
            testCase.verifyFalse(any([sosMeta.gram_basis_mismatch]));
        end

        function missingGramHandleFailsClosed(testCase)
            [~, ~, polyMeta, cfg] = SosPostcheckGramBindingRev3Test.buildMinimalProblem(4, 2);
            firstSos = find([polyMeta.sos_meta.is_sos_constraint], 1, 'first');
            polyMeta.sos_meta(firstSos).sos_constraint = [];

            oracle = independent_sos_postcheck_rev3(polyMeta, cfg, 0);
            testCase.verifyGreaterThanOrEqual(oracle.gram_handle_missing_count, 1);
            testCase.verifyFalse(oracle.strict_pass);
            testCase.verifyFalse(oracle.tol_pass);
        end

        function gramDimensionMismatchFailsClosed(testCase)
            [~, ~, polyMeta, cfg] = SosPostcheckGramBindingRev3Test.buildMinimalProblem(4, 2);
            firstSos = find([polyMeta.sos_meta.is_sos_constraint], 1, 'first');
            polyMeta.sos_meta(firstSos).gram_dim = polyMeta.sos_meta(firstSos).gram_dim + 1;

            oracle = independent_sos_postcheck_rev3(polyMeta, cfg, 0);
            testCase.verifyGreaterThanOrEqual(oracle.gram_dim_mismatch_count, 1);
            testCase.verifyFalse(oracle.strict_pass);
            testCase.verifyFalse(oracle.tol_pass);
        end
    end

    methods (Static)
        function [F, params, polyMeta, cfg] = buildMinimalProblem(dV, dmult)
            yalmip('clear');
            reg = SosPostcheckGramBindingRev3Test.boxRegion(-1, 1, -1, 1);
            Phi = {eye(2)};
            q = {[0; 0]};
            transitions = struct('i', {}, 'j', {}, 'G', {}, 'w', {});
            originRegions = [];
            cfg = SosPostcheckGramBindingRev3Test.cfg();
            degreeMeta = SosPostcheckGramBindingRev3Test.degreeMeta(dV, dmult);

            [F, params, polyMeta] = build_sos_problem_rev3(reg, Phi, q, transitions, ...
                originRegions, dV, dmult, degreeMeta, cfg);
        end

        function cfg = cfg()
            cfg = struct();
            cfg.epsilon = 1e-6;
            cfg.alpha = 1e-6;
            cfg.coefficient_bound = 1e3;
            cfg.normalization_mode = "fixed_margin_plus_box";
            cfg.gram_eig_tol = 1e-7;
            cfg.coeff_match_tol = 1e-2;
            cfg.sos_oracle_tol_factor = 10;
        end

        function ops = options()
            ops = sdpsettings('verbose', 0);
            if ~isempty(which('mosekopt'))
                ops.solver = 'mosek';
            end
        end

        function degreeMeta = degreeMeta(dV, dmult)
            residualDegree = max(dV, dmult + 1);
            dSigma = SosPostcheckGramBindingRev3Test.smallestEvenAtLeast(residualDegree);
            degreeMeta = struct();
            degreeMeta.d_V = dV;
            degreeMeta.d_mult = dmult;
            degreeMeta.residual_polynomial_degree = residualDegree;
            degreeMeta.d_sigma = dSigma;
            degreeMeta.residual_sos_degree = dSigma;
            degreeMeta.residual_has_odd_top_degree = mod(residualDegree, 2) ~= 0;
            degreeMeta.lyapunov_monomial_count = SosPostcheckGramBindingRev3Test.monomialCount2d(dV);
            degreeMeta.residual_gram_dim = SosPostcheckGramBindingRev3Test.monomialCount2d(dSigma / 2);
            if dmult == 0
                degreeMeta.multiplier_gram_dim = 1;
            else
                degreeMeta.multiplier_gram_dim = SosPostcheckGramBindingRev3Test.monomialCount2d(dmult / 2);
            end
        end

        function reg = boxRegion(xmin, xmax, ymin, ymax)
            reg = struct();
            reg.H = [1 0; -1 0; 0 1; 0 -1];
            reg.b = [xmax; -xmin; ymax; -ymin];
        end

        function d = smallestEvenAtLeast(v)
            d = ceil(v);
            if mod(d, 2) ~= 0
                d = d + 1;
            end
        end

        function n = monomialCount2d(d)
            n = nchoosek(2 + d, d);
        end
    end
end
