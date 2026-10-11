%% Run Rev3 SOS candidate monomial basis tests.

clear; clc;

test_dir = fileparts(mfilename('fullpath'));
release_root = find_release_root_for_sos_basis_tests(test_dir);
proposed_dir = fullfile(release_root, 'src', 'proposed');

addpath(proposed_dir);
addpath(fullfile(proposed_dir, 'functions'));
addpath(test_dir);
ensure_yalmip_paths_for_sos_basis_tests();

fprintf('Rev3 SOS basis tests\n');
fprintf('Proposed-method dir: %s\n', proposed_dir);

results = runtests('SosBasisRev3Test');
disp(results);
assert(all([results.Passed]), 'One or more Rev3 SOS basis tests failed.');


function ensure_yalmip_paths_for_sos_basis_tests()
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


function root = find_release_root_for_sos_basis_tests(start_dir)
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
