function D = build_time_delay_data()
%BUILD_TIME_DELAY_DATA Build the time-delay LMI representation from the
% canonical seven-region explicit PWA controller used in this package.

reg = answer_regions();
Ac = [0, 1; 0, 0];
Bc = [0; 1];
M = numel(reg);

D = struct();
D.description = "Double integrator with the canonical input-constrained explicit MPC law";
D.conversionMode = "euler-double-integrator";
D.M = M;
D.A = repmat({Ac}, M, 1);
D.a = repmat({[0; 0]}, M, 1);
D.B = Bc;
D.K = reshape({reg.K}, [], 1);
D.k = reshape({reg.k}, [], 1);
D.E = cell(M, 1);
D.e = cell(M, 1);

for i = 1:M
    % The LMI code uses E_i x + e_i >= 0.  The controller data use
    % H_i x <= b_i, hence E_i=-H_i and e_i=b_i.
    D.E{i} = -reg(i).H;
    D.e{i} = reg(i).b;
end
end
