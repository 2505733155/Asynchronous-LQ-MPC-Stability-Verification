%% Paper configuration: Rev1 PWQ Lyapunov SDP margin sweep
% Reproduces the complete PWQ comparison reported in the paper on the
% common sampled-PWA model. Both existing multiplier modes are retained.

setenv('PWQ_REV1_FAIR_MOSEK', '1');
setenv('PWQ_REV1_MOSEK_TOL', '1e-8');
setenv('PWQ_REV1_MARGIN_VALUES', ...
    '1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10');
setenv('PWQ_REV1_RESULT_FILENAME', 'pwq_paper_margin_sweep_results.csv');

run(fullfile(fileparts(mfilename('fullpath')), ...
    'run_region_piecewise_lyapunov_sdp_rev1.m'));
