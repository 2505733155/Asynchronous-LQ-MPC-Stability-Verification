%% Diagnostic wrapper: sampled-data PWA time-delay LMI route with MOSEK and state scaling
% Applies the coordinate change x = diag([20,5]) z before assembling the LMI
% diagnostics, while keeping the original physical explicit-MPC law.

setenv('TD_PWA_SOLVER', 'mosek');
setenv('TD_PWA_MOSEK_TOL', '1e-8');
setenv('TD_PWA_MOSEK_THREADS', '1');
setenv('TD_PWA_STATE_SCALE', '20,5');
setenv('TD_PWA_SUMMARY_CSV', ...
    'time_delay_paper_scaled20_5_results.csv');

run(fullfile(fileparts(mfilename('fullpath')), ...
    'compute_tau_max_di_mpc_strong_adj_scan.m'));
