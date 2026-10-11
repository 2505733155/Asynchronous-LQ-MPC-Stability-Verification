classdef SosBasisRev3Test < matlab.unittest.TestCase
    % Regression tests for explicit Rev3 SOS candidate monomial bases.

    methods (Test)
        function defaultDmult2UsesExpectedBases(testCase)
            polyMeta = SosBasisRev3Test.buildMinimalProblem(4, 2);
            dims = cellfun(@length, polyMeta.sos_candidate_monomials);
            meta = polyMeta.sos_meta;

            testCase.verifyEqual(dims, [3 3 3 3 6]);
            testCase.verifyEqual([meta.degree], [2 2 2 2 4]);
            testCase.verifyEqual([meta.gram_dim], [3 3 3 3 6]);
            testCase.verifyEqual(string({meta.kind}), ...
                ["multiplier" "multiplier" "multiplier" "multiplier" "residual"]);
            testCase.verifyFalse(polyMeta.degree_meta.residual_has_odd_top_degree);
        end

        function dmult4UsesExplicitSixthDegreeResidualBasis(testCase)
            polyMeta = SosBasisRev3Test.buildMinimalProblem(4, 4);
            dims = cellfun(@length, polyMeta.sos_candidate_monomials);
            meta = polyMeta.sos_meta;

            testCase.verifyEqual(dims, [6 6 6 6 10]);
            testCase.verifyEqual([meta.degree], [4 4 4 4 6]);
            testCase.verifyEqual([meta.gram_dim], [6 6 6 6 10]);
            testCase.verifyEqual(polyMeta.degree_meta.residual_polynomial_degree, 5);
            testCase.verifyEqual(polyMeta.degree_meta.residual_sos_degree, 6);
            testCase.verifyTrue(polyMeta.degree_meta.residual_has_odd_top_degree);
        end

        function candidateListMatchesSosConstraintCount(testCase)
            polyMeta = SosBasisRev3Test.buildMinimalProblem(4, 2);
            testCase.verifyEqual(numel(polyMeta.sos_candidate_monomials), ...
                nnz([polyMeta.sos_meta.is_sos_constraint]));
        end
    end

    methods (Static)
        function polyMeta = buildMinimalProblem(dV, dmult)
            yalmip('clear');
            reg = SosBasisRev3Test.boxRegion(-1, 1, -1, 1);
            Phi = {eye(2)};
            q = {[0; 0]};
            transitions = struct('i', {}, 'j', {}, 'G', {}, 'w', {});
            originRegions = [];
            cfg = SosBasisRev3Test.cfg();
            degreeMeta = SosBasisRev3Test.degreeMeta(dV, dmult);

            [~, ~, polyMeta] = build_sos_problem_rev3(reg, Phi, q, transitions, ...
                originRegions, dV, dmult, degreeMeta, cfg);
        end

        function cfg = cfg()
            cfg = struct();
            cfg.epsilon = 1e-6;
            cfg.alpha = 1e-6;
            cfg.coefficient_bound = 1e3;
            cfg.normalization_mode = "fixed_margin_plus_box";
        end

        function degreeMeta = degreeMeta(dV, dmult)
            residualDegree = max(dV, dmult + 1);
            dSigma = SosBasisRev3Test.smallestEvenAtLeast(residualDegree);
            degreeMeta = struct();
            degreeMeta.d_V = dV;
            degreeMeta.d_mult = dmult;
            degreeMeta.residual_polynomial_degree = residualDegree;
            degreeMeta.d_sigma = dSigma;
            degreeMeta.residual_sos_degree = dSigma;
            degreeMeta.residual_has_odd_top_degree = mod(residualDegree, 2) ~= 0;
            degreeMeta.lyapunov_monomial_count = SosBasisRev3Test.monomialCount2d(dV);
            degreeMeta.residual_gram_dim = SosBasisRev3Test.monomialCount2d(dSigma / 2);
            if dmult == 0
                degreeMeta.multiplier_gram_dim = 1;
            else
                degreeMeta.multiplier_gram_dim = SosBasisRev3Test.monomialCount2d(dmult / 2);
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
