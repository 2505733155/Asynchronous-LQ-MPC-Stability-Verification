%% KKT-lifted SOS global direct-attempt diagnostic with fixed state scaling [20,5]

setenv('KKT_SOS_GLOBAL_ATTEMPT', '1');
setenv('KKT_SOS_LIFT_SCALING', 'fixed_state_20_5');
setenv('KKT_SOS_OUTPUT_SUFFIX', '_scaled20_5');

run(fullfile(fileparts(mfilename('fullpath')), ...
    'run_kkt_lifted_sos_standard_mpc.m'));
