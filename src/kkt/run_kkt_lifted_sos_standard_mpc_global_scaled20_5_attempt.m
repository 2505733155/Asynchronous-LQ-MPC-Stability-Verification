%% Fixed [20,5] state scaling, preserving caller environment and output suffix.
run_fixed_state_scaling();

function run_fixed_state_scaling()
keys = {'KKT_SOS_GLOBAL_ATTEMPT','KKT_SOS_LIFT_SCALING','KKT_SOS_OUTPUT_SUFFIX'};
previous = cellfun(@getenv,keys,'UniformOutput',false);
cleanup = onCleanup(@() restore_environment(keys,previous)); %#ok<NASGU>
setenv(keys{1},'1');
setenv(keys{2},'fixed_state_20_5');
if isempty(previous{3}), setenv(keys{3},'_scaled20_5_qp4dp'); end
run(fullfile(fileparts(mfilename('fullpath')),'run_kkt_lifted_sos_standard_mpc.m'));
end

function restore_environment(keys,previous)
for i=1:numel(keys), setenv(keys{i},previous{i}); end
end
