%% Run Rev3 forward-invariance oracle tests.

clear; clc;

test_dir = fileparts(mfilename('fullpath'));
release_root = find_release_root_for_forward_invariance_tests(test_dir);
proposed_dir = fullfile(release_root, 'src', 'proposed');

addpath(proposed_dir);
addpath(fullfile(proposed_dir, 'functions'));
addpath(test_dir);
sample_data_dir = fullfile(release_root, 'data', 'double_integrator');
if isfolder(sample_data_dir)
    addpath(sample_data_dir);
end
ensure_matlab_linprog_precedence_for_forward_invariance_tests();

fprintf('Rev3 forward-invariance oracle tests\n');
fprintf('Proposed-method dir: %s\n', proposed_dir);
fprintf('linprog: %s\n', which('linprog'));

results = runtests('ForwardInvarianceOracleRev3Test');
disp(results);
assert(~any([results.Failed]), 'One or more Rev3 forward-invariance tests failed.');


function ensure_matlab_linprog_precedence_for_forward_invariance_tests()
optim_path = fullfile(matlabroot, 'toolbox', 'optim', 'optim');
if isfolder(optim_path)
    addpath(optim_path, '-begin');
end
end


function root = find_release_root_for_forward_invariance_tests(start_dir)
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
