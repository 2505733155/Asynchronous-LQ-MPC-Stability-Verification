classdef ForwardInvarianceOracleRev3Test < matlab.unittest.TestCase
    % Unit tests for verify_forward_invariance_union_rev3.

    methods (Test)
        function fullSpacePartitionBypassesLp(testCase)
            reg = ForwardInvarianceOracleRev3Test.boxRegion(-1, 1, -1, 1);
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({2 * eye(2)}, {[10; -3]});
            proof = verify_forward_invariance_union_rev3(reg, Phi, q, ...
                ForwardInvarianceOracleRev3Test.scaleCtx(true), ...
                ForwardInvarianceOracleRev3Test.cfg());

            testCase.verifyTrue(proof.verified);
            testCase.verifyFalse(proof.failed);
            testCase.verifyEqual(string(proof.status), "verified_full_state_space_partition");
            testCase.verifyEqual(string(proof.method), "full_state_space_partition");
            testCase.verifyEqual(proof.lp_count, 0);
        end

        function singleBoxIdentityIsInvariant(testCase)
            reg = ForwardInvarianceOracleRev3Test.boxRegion(-1, 1, -1, 1);
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({eye(2)}, {[0; 0]});
            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q);

            testCase.verifyFalse(proof.verified);
            testCase.verifyTrue(proof.tolerance_verified);
            testCase.verifyEqual(string(proof.status), "forward_invariance_verified_with_tolerance");
            testCase.verifyEqual(proof.lp_count, 4);
            testCase.verifyEqual(proof.max_escape_margin, 0, 'AbsTol', 1e-9);
        end

        function singleBoxContractionIsInvariant(testCase)
            reg = ForwardInvarianceOracleRev3Test.boxRegion(-1, 1, -1, 1);
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({0.5 * eye(2)}, {[0; 0]});
            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q);

            testCase.verifyFalse(proof.verified);
            testCase.verifyTrue(proof.tolerance_verified);
            testCase.verifyEqual(string(proof.status), "forward_invariance_verified_with_tolerance");
            testCase.verifyLessThanOrEqual(proof.max_escape_margin, -0.5 + 1e-9);
        end

        function singleBoxTranslationOutFails(testCase)
            reg = ForwardInvarianceOracleRev3Test.boxRegion(-1, 1, -1, 1);
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({eye(2)}, {[0.25; 0]});
            cfg = ForwardInvarianceOracleRev3Test.cfg();
            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q, cfg);

            testCase.verifyFalse(proof.verified);
            testCase.verifyTrue(proof.failed);
            testCase.verifyEqual(string(proof.status), "failed_forward_invariant_region_union");
            ForwardInvarianceOracleRev3Test.verifyCounterexample(testCase, reg, Phi, q, proof, cfg.forward_invariance_tol);
        end

        function nonconvexTwoBoxUnionIdentityIsInvariant(testCase)
            reg = ForwardInvarianceOracleRev3Test.disjointTwoBoxes();
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({eye(2), eye(2)}, {[0; 0], [0; 0]});
            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q);

            testCase.verifyFalse(proof.verified);
            testCase.verifyTrue(proof.tolerance_verified);
            testCase.verifyEqual(string(proof.status), "forward_invariance_verified_with_tolerance");
            testCase.verifyEqual(proof.lp_count, 32);
        end

        function nonconvexTwoBoxUnionGapJumpFails(testCase)
            reg = ForwardInvarianceOracleRev3Test.disjointTwoBoxes();
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({eye(2), eye(2)}, {[1.5; 0], [0; 0]});
            cfg = ForwardInvarianceOracleRev3Test.cfg();
            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q, cfg);

            testCase.verifyFalse(proof.verified);
            testCase.verifyTrue(proof.failed);
            testCase.verifyEqual(proof.counterexample_region, 1);
            ForwardInvarianceOracleRev3Test.verifyCounterexample(testCase, reg, Phi, q, proof, cfg.forward_invariance_tol);
        end

        function nonconvexTwoBoxUnionLegalRegionJumpPasses(testCase)
            reg = ForwardInvarianceOracleRev3Test.disjointTwoBoxes();
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({eye(2), eye(2)}, {[3; 0], [0; 0]});
            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q);

            testCase.verifyFalse(proof.verified);
            testCase.verifyTrue(proof.tolerance_verified);
            testCase.verifyEqual(string(proof.status), "forward_invariance_verified_with_tolerance");
        end

        function subToleranceEscapeIsNotStrictlyVerified(testCase)
            tau = 1e-6;
            reg = ForwardInvarianceOracleRev3Test.boxRegion(-1, 1, -1, 1);
            cfg = ForwardInvarianceOracleRev3Test.cfg(tau);

            [PhiSmall, qSmall] = ForwardInvarianceOracleRev3Test.maps({eye(2)}, {[0.5 * tau; 0]});
            proofSmall = ForwardInvarianceOracleRev3Test.runOracle(reg, PhiSmall, qSmall, cfg);
            testCase.verifyFalse(proofSmall.verified);
            testCase.verifyTrue(proofSmall.tolerance_verified);
            testCase.verifyEqual(string(proofSmall.status), "forward_invariance_verified_with_tolerance");
            testCase.verifyGreaterThan(proofSmall.max_escape_margin, 0);
            testCase.verifyLessThan(proofSmall.max_escape_margin, tau);

            [PhiLarge, qLarge] = ForwardInvarianceOracleRev3Test.maps({eye(2)}, {[20 * tau; 0]});
            proofLarge = ForwardInvarianceOracleRev3Test.runOracle(reg, PhiLarge, qLarge, cfg);
            testCase.verifyFalse(proofLarge.verified);
            testCase.verifyTrue(proofLarge.failed);
            ForwardInvarianceOracleRev3Test.verifyCounterexample(testCase, reg, PhiLarge, qLarge, proofLarge, tau);
        end

        function ambiguousLpContinuesThenReportsNotVerified(testCase)
            ForwardInvarianceOracleRev3Test.ambiguousThenActualLinprog('reset');
            reg = ForwardInvarianceOracleRev3Test.boxRegion(-1, 1, -1, 1);
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({eye(2)}, {[0; 0]});
            cfg = ForwardInvarianceOracleRev3Test.cfg();
            cfg.forward_invariance_linprog_fcn = @ForwardInvarianceOracleRev3Test.ambiguousThenActualLinprog;

            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q, cfg);

            testCase.verifyFalse(proof.verified);
            testCase.verifyFalse(proof.failed);
            testCase.verifyEqual(string(proof.status), "forward_invariance_not_verified");
            testCase.verifyEqual(proof.ambiguous_lp_count, 1);
            testCase.verifyGreaterThan(proof.lp_count, 1);
            testCase.verifyEqual(proof.counterexample_exitflag, 0);
        end

        function ambiguousLpDoesNotHideLaterCounterexample(testCase)
            ForwardInvarianceOracleRev3Test.ambiguousThenActualLinprog('reset');
            reg = ForwardInvarianceOracleRev3Test.boxRegion(-1, 1, -1, 1);
            [Phi, q] = ForwardInvarianceOracleRev3Test.maps({eye(2)}, {[0; 0.25]});
            cfg = ForwardInvarianceOracleRev3Test.cfg();
            cfg.forward_invariance_linprog_fcn = @ForwardInvarianceOracleRev3Test.ambiguousThenActualLinprog;

            proof = ForwardInvarianceOracleRev3Test.runOracle(reg, Phi, q, cfg);

            testCase.verifyFalse(proof.verified);
            testCase.verifyTrue(proof.failed);
            testCase.verifyEqual(string(proof.status), "failed_forward_invariant_region_union");
            testCase.verifyEqual(proof.ambiguous_lp_count, 1);
            testCase.verifyGreaterThan(proof.lp_count, 1);
            ForwardInvarianceOracleRev3Test.verifyCounterexample(testCase, reg, Phi, q, proof, cfg.forward_invariance_tol);
        end

        function realRev3SampleTimesRemainForwardInvariant(testCase)
            releaseRoot = ForwardInvarianceOracleRev3Test.releaseRoot();
            sampleDataDir = fullfile(releaseRoot, 'data', 'double_integrator');
            testCase.assumeTrue(exist(fullfile(sampleDataDir, 'answer_regions.m'), 'file') == 2, ...
                'Explicit-MPC sample data are not present in this checkout.');
            addpath(sampleDataDir);
            regOrig = answer_regions();
            reg = ForwardInvarianceOracleRev3Test.scaleRegions(regOrig, [20, 5]);
            cfg = ForwardInvarianceOracleRev3Test.cfg();
            scaleCtx = ForwardInvarianceOracleRev3Test.scaleCtx(false);

            sampleTimes = 0:0.1:1.5;
            for kk = 1:numel(sampleTimes)
                [Phi, q] = ForwardInvarianceOracleRev3Test.sampleMaps(regOrig, sampleTimes(kk), [20, 5]);
                proof = verify_forward_invariance_union_rev3(reg, Phi, q, scaleCtx, cfg);
                testCase.verifyFalse(proof.verified, sprintf('Ts=%g has no validated LP upper bound', sampleTimes(kk)));
                testCase.verifyTrue(proof.tolerance_verified, sprintf('Ts=%g should pass the tolerance check', sampleTimes(kk)));
                testCase.verifyEqual(string(proof.status), "forward_invariance_verified_with_tolerance");
            end
        end
    end

    methods (Static)
        function proof = runOracle(reg, Phi, q, cfg)
            if nargin < 4
                cfg = ForwardInvarianceOracleRev3Test.cfg();
            end
            proof = verify_forward_invariance_union_rev3(reg, Phi, q, ...
                ForwardInvarianceOracleRev3Test.scaleCtx(false), cfg);
        end

        function cfg = cfg(tau)
            if nargin < 1
                tau = 1e-6;
            end
            cfg = struct();
            cfg.forward_invariance_tol = tau;
        end

        function scaleCtx = scaleCtx(fullSpaceCovered)
            scaleCtx = struct();
            scaleCtx.partition_full_space_covered = fullSpaceCovered;
        end

        function [Phi, q] = maps(PhiIn, qIn)
            Phi = PhiIn;
            q = qIn;
        end

        function reg = boxRegion(xmin, xmax, ymin, ymax)
            reg = struct();
            reg.H = [1 0; -1 0; 0 1; 0 -1];
            reg.b = [xmax; -xmin; ymax; -ymin];
        end

        function reg = disjointTwoBoxes()
            reg(1) = ForwardInvarianceOracleRev3Test.boxRegion(-2, -1, -1, 1);
            reg(2) = ForwardInvarianceOracleRev3Test.boxRegion(1, 2, -1, 1);
        end

        function verifyCounterexample(testCase, reg, Phi, q, proof, tau)
            i = proof.counterexample_region;
            testCase.verifyGreaterThanOrEqual(i, 1);
            testCase.verifyLessThanOrEqual(i, numel(reg));
            x = proof.counterexample_x;
            xp = proof.counterexample_xp;

            testCase.verifyLessThanOrEqual(max(reg(i).H * x - reg(i).b), tau);
            testCase.verifyEqual(xp, Phi{i} * x + q{i}, 'AbsTol', 1e-9);

            for j = 1:numel(reg)
                escapeMargin = max(reg(j).H * xp - reg(j).b);
                testCase.verifyGreaterThanOrEqual(escapeMargin, 0.5 * tau);
            end
        end

        function root = releaseRoot()
            thisFile = mfilename('fullpath');
            testDir = fileparts(thisFile);
            root = testDir;
            while true
                if exist(fullfile(root, 'data', 'double_integrator', 'answer_regions.m'), 'file') == 2
                    return;
                end
                parent = fileparts(root);
                if strcmp(parent, root)
                    error('Could not find release root from %s', testDir);
                end
                root = parent;
            end
        end

        function reg = scaleRegions(regOrig, stateScale)
            T = diag(stateScale);
            reg = regOrig;
            for i = 1:numel(regOrig)
                reg(i).H = regOrig(i).H * T;
                reg(i).b = regOrig(i).b;
            end
        end

        function [Phi, q] = sampleMaps(regOrig, Ts, stateScale)
            M = numel(regOrig);
            Phi = cell(1, M);
            q = cell(1, M);
            A0 = [1 Ts; 0 1];
            B0 = [Ts^2 / 2; Ts];
            T = diag(stateScale);

            for i = 1:M
                PhiOriginal = A0 + B0 * regOrig(i).K;
                qOriginal = B0 * regOrig(i).k;
                Phi{i} = T \ (PhiOriginal * T);
                q{i} = T \ qOriginal;
            end
        end

        function [x, fval, exitflag] = ambiguousThenActualLinprog(varargin)
            persistent callCount
            if nargin == 1 && (ischar(varargin{1}) || isstring(varargin{1})) && string(varargin{1}) == "reset"
                callCount = 0;
                x = [];
                fval = [];
                exitflag = [];
                return;
            end
            if isempty(callCount)
                callCount = 0;
            end
            callCount = callCount + 1;
            if callCount == 1
                f = varargin{1};
                x = NaN(size(f));
                fval = NaN;
                exitflag = 0;
                return;
            end
            [x, fval, exitflag] = linprog(varargin{:});
        end
    end
end
