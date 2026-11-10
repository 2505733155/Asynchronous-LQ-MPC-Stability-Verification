function qp = double_integrator_qp_data()
%DOUBLE_INTEGRATOR_QP_DATA Numerical QP underlying the archived PWA example.
% Objective: 0.5*U'*H*U + U'*Fux*x, with -1 <= U_i <= 1.
% Fxu is the transpose used in the manuscript's x'*F*U convention.
% Nominal design metadata must not be used to overwrite these QP matrices.
qp = struct();
qp.H = [0.8365 0.3603;0.3603 0.2059];
qp.Fux = [0.4624 1.2852;0.1682 0.5285];
qp.Fxu = qp.Fux';
qp.N = 2;
qp.u_min = -1;
qp.u_max = 1;
qp.source_id = "hybtbx_normalized_qp_4dp";
qp.nominal = struct('Ap',[1 1;0 1],'Bp',[0;1], ...
    'Tsp',1,'Q',diag([1 0]),'R',0.1,'N',2, ...
    'terminal_weight',"DARE for the nominal discrete model");
end
