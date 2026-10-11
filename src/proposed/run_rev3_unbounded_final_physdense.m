%% Rev3 unbounded input-only explicit-MPC SOS scan for the selected setting.
% Selected diagnostic setting:
%   gamma = epsilon = alpha = 1e-2
%   MOSEK primal/dual/gap tolerance = 1e-8
%   dense postcheck physical box: x1 in [-100,100], x2 in [-100,100]

setenv('REV3_RUN_LABEL', 'unbounded_gamma1em2_tol1em8_physdense_x100_v100_15x57');
setenv('REV3_OUTPUT_STEM', 'region_piecewise_polynomial_sos_rev3_unbounded_gamma1em2_tol1em8_physdense_x100_v100');
setenv('REV3_SAMPLE_TIMES', '0.1:0.1:1.0');
setenv('REV3_EPSILON', '1e-2');
setenv('REV3_ALPHA', '1e-2');
setenv('REV3_MOSEK_TOL', '1e-8');
setenv('REV3_DENSE_PHYSICAL_BOX', '-100 100 -100 100');
setenv('REV3_DENSE_GRID_COUNTS', '15 57');

run(fullfile(fileparts(mfilename('fullpath')), ...
    'run_region_piecewise_polynomial_sos_rev3.m'));
