function test_kkt_wrapper_environment()
%TEST_KKT_WRAPPER_ENVIRONMENT Exercise the real wrapper on normal/error exits.
root=fileparts(fileparts(mfilename('fullpath')));
source=fullfile(root,'src','kkt','run_kkt_lifted_sos_standard_mpc_global_scaled20_5_attempt.m');
keys={'KKT_SOS_GLOBAL_ATTEMPT','KKT_SOS_LIFT_SCALING','KKT_SOS_OUTPUT_SUFFIX'};
original=cellfun(@getenv,keys,'UniformOutput',false);
tempRoot=tempdir; target=tempname(tempRoot); mkdir(target);
cleanup=onCleanup(@() restore_and_remove(keys,original,target,tempRoot)); %#ok<NASGU>
copyfile(source,fullfile(target,'run_kkt_lifted_sos_standard_mpc_global_scaled20_5_attempt.m'));
for i=1:3, setenv(keys{i},['sentinel_' num2str(i)]); end
before=cellfun(@getenv,keys,'UniformOutput',false);
stub=fullfile(target,'run_kkt_lifted_sos_standard_mpc.m');
fid=fopen(stub,'w'); assert(fid>=0);
fprintf(fid,"assert(strcmp(getenv('KKT_SOS_GLOBAL_ATTEMPT'),'1'));\n");
fprintf(fid,"assert(strcmp(getenv('KKT_SOS_LIFT_SCALING'),'fixed_state_20_5'));\n");
fprintf(fid,"assert(strcmp(getenv('KKT_SOS_OUTPUT_SUFFIX'),'sentinel_3'));\n");
fprintf(fid,"setenv('KKT_SOS_OUTPUT_SUFFIX','changed_by_stub');\n");
fclose(fid);
run(fullfile(target,'run_kkt_lifted_sos_standard_mpc_global_scaled20_5_attempt.m'));
assert(isequal(before,cellfun(@getenv,keys,'UniformOutput',false)));
fid=fopen(stub,'w'); assert(fid>=0);
fprintf(fid,"setenv('KKT_SOS_OUTPUT_SUFFIX','changed_before_error');\nerror('test:ExpectedWrapperError','controlled error');\n");
fclose(fid);
try
    run(fullfile(target,'run_kkt_lifted_sos_standard_mpc_global_scaled20_5_attempt.m'));
    error('test:MissingExpectedError','The controlled failure did not occur.');
catch err
    assert(strcmp(err.identifier,'test:ExpectedWrapperError'));
end
assert(isequal(before,cellfun(@getenv,keys,'UniformOutput',false)));
fprintf('PASS: original wrapper preserves explicit suffix and restores environment on success/error.\n');
end

function restore_and_remove(keys,original,target,tempRoot)
for i=1:numel(keys),setenv(keys{i},original{i});end
% Delete only this test's freshly created directory within the system temp root.
assert(startsWith(target,tempRoot) && ~strcmp(target,tempRoot));
if isfolder(target),rmdir(target,'s');end
end
