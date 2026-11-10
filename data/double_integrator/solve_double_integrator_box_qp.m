function [U, ok, residual] = solve_double_integrator_box_qp(x, qp)
%SOLVE_DOUBLE_INTEGRATOR_BOX_QP Exhaust all nine strictly convex box-QP modes.
% Independent of quadprog and its possible MOSEK path shadowing.
H=qp.H; f=qp.Fux*x(:);
assert(isequal(size(H),[2 2]) && numel(f)==2 && min(eig(H))>0);
assert(qp.u_min==-1 && qp.u_max==1);
U=NaN(2,1); residual=Inf; best=Inf; tolerance=1e-10;
for first=-1:1
    for second=-1:1
        mode=[first;second]; free=(mode==0); candidate=mode;
        candidate(free)=-H(free,free)\(f(free)+H(free,~free)*candidate(~free));
        gradient=H*candidate+f;
        rv=max([0;abs(candidate)-1;abs(gradient(free)); ...
            -gradient(mode==-1);gradient(mode==1)]);
        if rv>tolerance, continue; end
        value=0.5*candidate'*H*candidate+f'*candidate;
        if value<best
            best=value; U=candidate; residual=rv;
        end
    end
end
ok=isfinite(best) && all(isfinite(U)) && residual<=tolerance;
end
