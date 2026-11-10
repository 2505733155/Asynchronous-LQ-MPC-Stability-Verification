function test_double_integrator_qp_data()
%TEST_DOUBLE_INTEGRATOR_QP_DATA Check the numerical QP independently of SDP.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'data','double_integrator'));
assert(qp4dp_result_suffix('')=="_qp4dp");
assert(qp4dp_result_suffix('_scaled20_5')=="_scaled20_5_qp4dp");
assert(qp4dp_result_suffix('_small_box')=="_small_box_qp4dp");
assert(qp4dp_result_suffix('_scaled20_5_qp4dp_m04')=="_scaled20_5_qp4dp_m04");
qp = double_integrator_qp_data();
assert(isequal(qp.H,[0.8365 0.3603;0.3603 0.2059]));
assert(isequal(qp.Fux,[0.4624 1.2852;0.1682 0.5285]));
assert(isequal(qp.Fxu,qp.Fux') && min(eig(qp.H))>0);
K = -qp.H\qp.Fux;
assert(norm(K(1,:)-[-0.815801595784556 -1.749302793118031])<2e-14);
assert(norm(K(1,:)-[-0.816617139036083 -1.749930783859335])>1e-4);
[U,ok,residual] = solve_double_integrator_box_qp([0.1;0],qp);
assert(ok && residual<1e-12 && norm(U-K*[0.1;0])<1e-13);
reg = answer_regions_input_only_legacy();
rng(109,'twister');
points = [zeros(2,1),[50;-20],[-50;20],40*rand(2,300)-20];
for x=points
    [U,ok,residual] = solve_double_integrator_box_qp(x,qp);
    assert(ok && residual<1e-10 && max(abs(U))<=1+1e-12);
    grad=qp.H*U+qp.Fux*x;
    assert(all(grad(U<-1+1e-10)>=-1e-10));
    assert(all(grad(U>1-1e-10)<=1e-10));
    found=false;
    for i=1:numel(reg)
        if all(reg(i).H*x<=reg(i).b+1e-10)
            assert(abs(reg(i).K*x+reg(i).k-U(1))<2e-12);
            found=true;
        end
    end
    assert(found);
end
fprintf('PASS: QP data, F orientation, 303 state/KKT checks and legacy inputs.\n');
end
