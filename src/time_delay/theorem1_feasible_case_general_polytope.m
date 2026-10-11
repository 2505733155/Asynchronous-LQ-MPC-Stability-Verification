function [ok, cert] = theorem1_feasible_case_general_polytope(tau, A, a, B, K, k, E, e, I0, nx, opts, enforceEtaGamma, enforceLambdaEntrywise, etaGammaScale, gammaLowerBound)
% Feasibility check for Theorem 1 LMIs at fixed tau.
% Generalized from example/theorem1_feasible_case.m to support arbitrary
% numbers of region inequalities p_i (not restricted to p_i = nx).

epsPD = 1e-7;
M = numel(A);
nu = size(B,2);

if nargin < 12 || isempty(enforceEtaGamma)
    enforceEtaGamma = false;
end
if nargin < 13 || isempty(enforceLambdaEntrywise)
    enforceLambdaEntrywise = true;
end
if nargin < 14 || isempty(etaGammaScale)
    etaGammaScale = 1;
end
if nargin < 15 || isempty(gammaLowerBound)
    gammaLowerBound = 1 + 1e-6;
end

if isempty(I0)
    notI0 = 1:M;
else
    notI0 = setdiff(1:M, I0);
end

yalmip('clear');

P = sdpvar(nx, nx, 'symmetric');
R = sdpvar(nx, nx, 'symmetric');
X = sdpvar(nx, nx, 'symmetric');
eta = sdpvar(1,1);
gamma = sdpvar(1,1);

N = cell(M,1);
for i = 1:M
    N{i} = sdpvar(2*nx, nx, 'full');
end

Lambda = cell(M,1);
for i = notI0
    pi_i = size(E{i},1);
    Lambda{i} = sdpvar(pi_i + 1, pi_i + 1, 'symmetric');
end

C = [];
C = [C, P >= epsPD*eye(nx), R >= epsPD*eye(nx), X >= epsPD*eye(nx)];
if enforceEtaGamma
    C = [C, gamma >= gammaLowerBound, eta >= etaGammaScale * gamma];
else
    C = [C, eta >= epsPD, gamma >= epsPD];
end

Icol = [eye(nx); -eye(nx)];
Irow = [eye(nx) -eye(nx)];

for i = 1:M
    Ai = A{i};
    Ki = K{i};
    BKi = B * Ki;

    if ~isempty(I0) && any(i == I0)
        Ni = N{i};

        Psi = [P; zeros(nx,nx)]*[Ai BKi] + [Ai'; BKi']*[P zeros(nx,nx)] ...
            - Icol*X*Irow - Ni*Irow - Icol*Ni' + eta*eye(2*nx);

        Omega = [Psi, [P; zeros(nx,nx)]*B; B'*[P zeros(nx,nx)], -gamma*eye(nu)];

        M1 = [Ai'; BKi'; B'] * R * [Ai BKi B];
        M2 = [Icol; zeros(nu,nx)] * X * [Ai BKi B] ...
            + [Ai'; BKi'; B'] * X * [Irow zeros(nx,nu)];

        L18 = Omega + tau*M1 + tau*M2;
        C = [C, L18 <= -epsPD*eye(size(L18,1))];

        Q19 = [Omega + tau*M1, tau*[Ni; zeros(nu,nx)]; ...
            tau*[Ni' zeros(nx,nu)], -tau*R];
        C = [C, Q19 <= -epsPD*eye(size(Q19,1))];

    else
        Ni = N{i};
        Li = Lambda{i};
        pi_i = size(E{i},1);
        ci = B*k{i} + a{i};

        Ibar_col = [Icol; zeros(1,nx)];
        Ibar_row = [Irow zeros(nx,1)];
        Nbar_col = [Ni; zeros(1,nx)];
        Nbar_row = [Ni' zeros(nx,1)];

        Sleft  = [E{i}' zeros(nx,1); zeros(nx,pi_i) zeros(nx,1); e{i}' 1];
        % General-polytope fix: second row block uses zeros(1,nx), not zeros(1,pi_i).
        Sright = [E{i} zeros(pi_i,nx) e{i}; zeros(1,nx) zeros(1,nx) 1];

        barPsi = [P; zeros(nx,nx); zeros(1,nx)]*[Ai BKi ci] ...
            + [Ai'; BKi'; ci']*[P zeros(nx,nx) zeros(nx,1)] ...
            - Ibar_col*X*Ibar_row - Nbar_col*Ibar_row - Ibar_col*Nbar_row ...
            + blkdiag(eta*eye(2*nx), 0) + Sleft*Li*Sright;

        barOmega = [barPsi, [P; zeros(nx,nx); zeros(1,nx)]*B; ...
            B'*[P zeros(nx,nx) zeros(nx,1)], -gamma*eye(nu)];

        barM1 = [Ai'; BKi'; ci'; B'] * R * [Ai BKi ci B];
        barM2 = [Icol; zeros(1,nx); zeros(nu,nx)] * X * [Ai BKi ci B] ...
            + [Ai'; BKi'; ci'; B'] * X * [Irow zeros(nx,1) zeros(nx,nu)];

        L20 = barOmega + tau*barM1 + tau*barM2;
        C = [C, L20 <= -epsPD*eye(size(L20,1))];

        Q21 = [barOmega + tau*barM1, tau*[Ni; zeros(1 + nu,nx)]; ...
            tau*[Ni' zeros(nx,1 + nu)], -tau*R];
        C = [C, Q21 <= -epsPD*eye(size(Q21,1))];

        if enforceLambdaEntrywise
            C = [C, Li >= epsPD*eye(pi_i + 1), Li(:) >= 0];
        else
            C = [C, Li >= epsPD*eye(pi_i + 1)];
        end
    end
end

sol = optimize(C, [], opts);
ok = (sol.problem == 0);

if nargout > 1
    cert = struct();
    cert.problem = sol.problem;
    cert.info = sol.info;
    cert.tau = tau;
    cert.etaGammaScale = etaGammaScale;
    cert.gammaLowerBound = gammaLowerBound;
    cert.P = [];
    cert.R = [];
    cert.X = [];
    cert.eta = [];
    cert.gamma = [];
    cert.N = [];
    cert.Lambda = [];

    if ok
        cert.P = value(P);
        cert.R = value(R);
        cert.X = value(X);
        cert.eta = value(eta);
        cert.gamma = value(gamma);

        cert.N = cell(M,1);
        for i = 1:M
            cert.N{i} = value(N{i});
        end

        cert.Lambda = cell(M,1);
        for i = notI0
            cert.Lambda{i} = value(Lambda{i});
        end
    end
end
end
