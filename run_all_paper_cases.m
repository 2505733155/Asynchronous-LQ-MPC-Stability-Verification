function run_all_paper_cases()
%RUN_ALL_PAPER_CASES Run the four paper configurations.
% This sequential run can take substantial time and requires the external
% dependencies listed in README.md.

release_root = fileparts(mfilename('fullpath'));
data_dir = fullfile(release_root, 'data', 'double_integrator');
addpath(data_dir);
validate_double_integrator_data();

runners = {
    fullfile('src', 'proposed', 'run_rev3_unbounded_final_physdense.m')
    fullfile('src', 'pwq', 'run_pwq_rev1_mosek_paper_margin_sweep.m')
    fullfile('src', 'time_delay', 'run_time_delay_mosek_scaled_diag.m')
    fullfile('src', 'kkt', 'run_kkt_lifted_sos_standard_mpc_global_scaled20_5_attempt.m')
    };

for i = 1:numel(runners)
    fprintf('\n===== Running paper configuration %d of %d: %s =====\n', ...
        i, numel(runners), runners{i});
    run_route(release_root, runners{i});
end

fprintf('\nAll paper configurations completed. Generated files are in results/generated/.\n');
end


function run_route(release_root, relative_runner)
runner = fullfile(release_root, relative_runner);
if exist(runner, 'file') ~= 2
    error('Missing paper-configuration runner: %s', runner);
end
% Each route is intentionally a script and begins by clearing its own
% workspace. Isolating it in this helper keeps the orchestrator state.
run(runner);
end
