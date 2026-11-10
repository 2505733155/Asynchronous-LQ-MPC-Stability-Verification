"""Portable rational verification of the HMPC original-proof comparison.

Reads the package actual QP/PWA data and retains both actual gains through an
explicit input-history feedback bound. Uses Python standard library only.
The fixed-reference proof was independently checked; the reported ranges are
sufficient inner estimates, not maximum capabilities of HMPC.
"""
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json
import re
import runpy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PUBLIC = ROOT
HELPER = PUBLIC / 'src/hmpc/verify_hmpc_reference.py'
A = runpy.run_path(str(HELPER))  # Arithmetic functions only; do not call calculate/W(h).
mm, tr, add, scale = (A[k] for k in ('matmul', 'transpose', 'add', 'scale'))
inv, pos, eye = (A[k] for k in ('inverse2', 'positive2', 'eye'))
PROFILES = {
    'higher_disturbance': ('.11', '.85', '.1', '.02', '22', '.0005'),
    'balanced': ('.12', '.94', '.04', '.01', '48', '.00025'),
    'higher_period': ('.13', '.99', '.005', '.002', '377', '.00003'),
}


def sqrt_bounds(value):
    assert value >= 0
    grid = 10**12
    n = isqrt(value.numerator * grid**2 // value.denominator)
    lo = F(n, grid)
    hi = lo if lo*lo == value else F(n+1, grid)
    assert lo*lo <= value <= hi*hi
    return lo, hi


def det(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def eigen_bounds(a, metric=None):
    # Eigenvalues of a symmetric pair a v = lambda metric v.
    metric = eye(2) if metric is None else metric
    assert pos(metric)
    transformed = mm(inv(metric), a)
    trace = transformed[0][0]+transformed[1][1]
    discriminant = trace*trace-4*det(a)/det(metric)
    root_hi = sqrt_bounds(discriminant)[1]
    return (trace-root_hi)/2, (trace+root_hi)/2


def row_square(row, metric_inverse):
    return mm([row], mm(metric_inverse, tr([row])))[0][0]


def read_fixed_data():
    pc = [[F('2.1429023470226376'), F('1.2245640609262458')],
          [F('1.2245640609262458'), F('1.3995571393121802')]]
    qc, rc = [[F(1), F(0)], [F(0), F('.0001')]], F('.1')
    qp = PUBLIC / 'data/double_integrator/double_integrator_qp_data.m'
    law = PUBLIC / 'data/double_integrator/answer_regions_input_only_legacy.m'
    text, lawtext = qp.read_text(encoding='utf-8'), law.read_text(encoding='utf-8')

    def matrix(name):
        v = re.search(r'qp\.'+name+r'\s*=\s*\[([^]]+)\]', text)[1]
        return [[F(x) for x in row.split()] for row in v.split(';')]

    lu = scale(mm(inv(matrix('H')), matrix('Fux')), -1)
    stored = [F(x.strip()) for x in re.search(r'K\s*=\s*\{\s*\[([^]]+)\]', lawtext)[1].split(',')]
    origin = re.search(r'H\s*=\s*\{\s*\[([^]]+)\]', lawtext)[1]
    origin = [[F(x.strip()) for x in row.split(',')] for row in origin.split(';')]
    kd = A['discrete_gain'](pc, qc, rc, 2)
    return pc, qc, rc, kd, lu, stored, origin, qp, law


def reference_interface(pc, qc, rc, kd, actuals):
    ct, ctbox, _ = A['reference_ct_enclosure'](pc, qc)
    stored = actuals['stored_PWA']
    e = [stored[0]-F('-3.1973071'), stored[1]-F('-2.5482575')]
    metric = add(mm(tr([e]), [e]), scale(eye(2), F('.0001')))
    im = inv(metric)
    for corner in product(*ctbox):
        assert row_square([x-y for x,y in zip(kd, corner)], im) < F('1.01')**2
    for actual in actuals.values():
        assert row_square([x-y for x,y in zip(actual,kd)], im) < F('.04')**2
    # Recheck the previously independently reviewed CT convolution estimate.
    assert F('1.274') < -ct[1][1]/2 <= -ct[1][0]/2 < F('1.275')
    assert F('1.25')**2 < -ct[0][1]-F('1.275')**2
    assert -ct[0][0]-F('1.274')**2 < F('1.26')**2
    bhi = (e[0]-F('1.27')*e[1])/F('1.25')
    blo = (e[0]-F('1.28')*e[1])/F('1.26')
    assert blo > e[1] > 0 and e[1]**2+bhi**2 < F('1.37')**2
    assert 3*F('3.14')/(4*F('1.26')) > F('1.86')
    expb = A['exp_positive_bounds']
    assert expb(F('1.27')*F('1.86'))[0] > 1/F('.095')
    assert expb(F('1.27')*F('3.14')/F('1.26'))[0] > 1/F('.043')
    assert (1/F('1.25'))**2+(1+F('1.28')/F('1.25'))**2 < F('2.2')**2
    gamma = e[0]/F('3.19')+2*F('1.37')*F('1.26')*F('.095')/(F('3.19')*(1-F('.043')))+F('.01')*F('2.2')/F('1.27')
    assert gamma < F('.88')
    assert F('.88')*F('1.01') == F('.8888') < 1
    gamma_d = F('.88')/(1-F('.8888'))
    assert gamma_d*F('.04') < 1
    return {'CT_gain_upper': .88, 'CT_to_Euler_error_upper': 1.01,
            'CT_to_Euler_product': .8888,
            'reference_to_actual_continuous_error_product': float(gamma_d*F('.04')),
            'norm_scope': 'These CT-stage constants use S, not the later Y metric.',
            'family_scope': 'Same Pc/Omega and genuine Euler family d=2/N. Rounded QP is NOT asserted identical to it.'}


def calculate(period=None, level=F('.00085'), metric=None, profile='higher_period'):
    hp,rho_theta,rho_input,eps,gamma_reference,disturbance = map(F,PROFILES[profile])
    period = hp if period is None else period
    assert 0 < period <= F(1) and level > 0
    y = [[F(1),F('.98')],[F('.98'),F('1.4272')]] if metric is None else metric
    assert pos(y), 'ISS metric must be positive definite'
    iy = inv(y)
    pc, qc, rc, kd, lu, stored, origin, qp, law = read_fixed_data()
    actuals = {'ideal_rounded_QP':lu[0], 'stored_PWA':stored}
    interface = reference_interface(pc,qc,rc,kd,actuals)
    ac, bc = [[F(0),F(1)],[F(0),F(0)]], [[F(0)],[F(1)]]
    kf = scale(mm(tr(bc),pc),-1/rc)
    terminal_residual = add(qc,add(mm(tr(ac),pc),add(mm(pc,ac),scale(mm(mm(pc,bc),mm(tr(bc),pc)),-1/rc))))
    assert pos(scale(terminal_residual,-1)) and F('.003')*pc[1][1]/rc**2 < 1

    # Sampling stage: original (15)-(18), written in the Y norm, with G=1.
    ak = add(ac,mm(bc,[kd]))
    dy = scale(add(mm(tr(ak),y),mm(y,ak)),-1)
    q, _ = eigen_bounds(dy,y)
    assert q > 0 and pos(add(dy,scale(y,-q)))
    pb = sqrt_bounds(y[1][1])[1]
    a,b = -kd[0],-kd[1]
    cx = sqrt_bounds(row_square([F(0),-a],iy))[1]+b*sqrt_bounds(row_square(kd,iy))[1]
    ctheta,cinput = 2*pb/(q*rho_theta),2*pb/(q*rho_input)
    alpha3 = 2*q*(1-rho_theta-rho_input)
    assert alpha3 > 0
    assert q-2*pb/cinput >= alpha3
    assert q-2*pb/ctheta-2*pb/cinput == alpha3/2
    gamma1,gamma2 = (1+eps)*cx,(1+1/eps)*b
    multiplier = ctheta*period*gamma1
    assert multiplier < 1, 'This metric/gain witness fails original equation (18)'
    assert gamma_reference >= max(ctheta*period*gamma2,cinput)
    error_upper = F('.0008')
    errors = {}
    for name,k in actuals.items():
        error = [x-v for x,v in zip(k,kd)]
        assert row_square(error,iy) < error_upper**2
        errors[name] = float(sqrt_bounds(row_square(error,iy))[1])
    actual_feedback_gain = 2*gamma_reference*error_upper
    assert actual_feedback_gain < 1
    # Additive-input interface: delta3=e(x[k])+w. Its history extends to 3h;
    # apply the SAME original Razumikhin Theorem 1 with max gain below one.
    combined_history_gain = max(multiplier,actual_feedback_gain)
    assert combined_history_gain < 1

    # Domain stage: bounding policy only. Pbar is NOT a new terminal cost.
    pbar = scale(pc,F('1.01'))
    abar = add(ac,mm(bc,[stored]))
    md = add(eye(2),abar)
    stage = add(qc,scale(mm(tr([stored]),[stored]),rc))
    delta1 = add(pbar,scale(add(mm(tr(md),mm(pbar,md)),stage),-1))
    assert pos(delta1) and pos(add(pbar,scale(pc,-1)))
    # At any 0<d<=1 the Euler residual is
    # d*[delta1+(1-d)*abar' Pbar abar] > 0, so ALL N>=2 are covered.
    assert F('.003')*pbar[1][1]/rc**2 < 1
    rho,mu = F('.003'),F('3.522')
    assert pos(add(scale(y,mu),scale(pbar,-1)))
    assert level*mu < rho, 'Level exceeds the proved common CT/DT valid domain'
    for row in lu+origin:
        assert level*row_square(row,iy) < 1
    eig_lo,eig_hi = eigen_bounds(y)
    assert eig_lo > 0
    chi = F('.0195')
    assert chi*chi*eig_hi < level
    # Original Theorem 1 max estimate, with the profile's 2*gamma input gain.
    input_gain = 2*gamma_reference
    assert input_gain**2*disturbance**2 < level
    for k in actuals.values():
        assert sqrt_bounds(level*row_square(k,iy))[1]+disturbance < 1
    cases=[]
    for j in range(1,11):
        h=F(j,10)
        accepted=(h<=period and ctheta*h*gamma1<1 and actual_feedback_gain<1)
        cases.append({'h':float(h),'this_witness_eq18_multiplier':float(ctheta*h*gamma1),
                      'candidate_original_route_pass':accepted,
                      'initial_Euclidean_ball':float(chi) if accepted else None})
    sources=[Path(__file__).resolve(),HELPER,qp,law]
    return {'version':'hmpc-original-proof-2026-11-10','profile':profile,
            'status':'Original-proof instantiation and actual-feedback closure independently verified; see technical derivation.',
            'literal_original_theorem_for_actual_QP':False,
            'classification':'Quantitative HMPC/Teel proof instantiation with actual-feedback input-history closure; not a literal equality of OCPs.',
            'no_W_equation_used':True,'no_new_sampling_theorem':True,
            'reference_interface':interface,
            'fixed_reference':{'T_H':2,'N':2,'d':1,'Qc':[[1,0],[0,.0001]],'Rc':.1,
                               'Pc':[[str(v) for v in row] for row in pc],'Omega_level':.003},
            'reference_gain_exact':[str(v) for v in kd],
            'actual_gain_exact':{k:[str(v) for v in row] for k,row in actuals.items()},
            'sampling':{'Y':[[str(v) for v in row] for row in y],
                        'q_lower':str(q),'PB_Y_norm_upper':str(pb),'history_slope_C_upper':str(cx),
                        'theta_allocation':str(rho_theta),'input_allocation':str(rho_input),'split_epsilon':str(eps),
                        'alpha3_coefficient':str(alpha3),'gamma_theta_slope':str(ctheta),'gamma_input_slope':str(cinput),
                        'gamma1_slope':str(gamma1),'gamma2_slope':str(gamma2),
                        'this_witness_strict_period_threshold':float(1/(ctheta*gamma1)),
                        'checked_fixed_period_interval_upper':float(period),
                        'eq18_multiplier_at_upper':float(multiplier),
                        'reference_ISS_gain_Y_upper':int(gamma_reference),
                        'actual_reference_error_Y_upper':float(error_upper),
                        'actual_reference_errors_Y_diagnostic':errors,
                        'actual_feedback_history_gain_upper':float(actual_feedback_gain),
                        'combined_history_gain_upper':float(combined_history_gain),'history_length':'3h after actual-input feedback closure'},
            'sampled_actual_feedback_error_closed':True,
            'domain':{'bounding_matrix_Pbar':'1.01 Pc, proof bound only; terminal Pc and Omega unchanged',
                      'Euler_d1_policy_residual':[[str(v) for v in row] for row in delta1],
                      'uniform_CT_and_all_Euler_N_ge_2_domain':'x^T Pbar x < .003',
                      'prediction_input_square_upper':float(rho*pbar[1][1]/rc**2),
                      'Pbar_le_mu_Y':float(mu),'Y_level':float(level),
                      'Y_eigenvalue_upper':str(eig_hi),
                      'initial_Euclidean_ball':float(chi),'input_disturbance_bound':float(disturbance),
                      'physical_trajectory_gain_upper':float(input_gain/sqrt_bounds(eig_lo)[0]),
                      'zero_disturbance_statement':'For fixed h in the checked interval: x0^T Y x0<=level remains in it and x(t)->0; no within-hold monotonicity is claimed.',
                      'scope_upper_bounds':{'fixed_original_terminal_set':.720081,'fixed_cost_horizon_any_admissible_terminal_design':2.08},
                      'upper_bound_provenance':'Previously derived necessary reference-feasibility bounds, not maximum actual ROA or all-HMPC bounds.'},
            'original_route_candidate_grid_points':[v['h'] for v in cases if v['candidate_original_route_pass']],
            'ten_points':cases,'controller_data_identity_checked':True,
            'sources':[{'path':str(v.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for v in sources]}



if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',choices=list(PROFILES),default='higher_period')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    result=calculate(profile=args.profile)
    if args.check:
        assert result==json.loads(args.check.read_text(encoding='utf-8')), 'Saved result or source identity differs'
    if args.output:
        assert not args.output.exists(), 'Do not overwrite existing evidence'
        args.output.write_bytes((json.dumps(result,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
    print(json.dumps({'profile':result['profile'],
                      'fixed_period_upper':result['sampling']['checked_fixed_period_interval_upper'],
                      'initial_Euclidean_ball':result['domain']['initial_Euclidean_ball'],
                      'input_disturbance_bound':result['domain']['input_disturbance_bound'],
                      'tested_grid_points_covered':result['original_route_candidate_grid_points']},indent=2))
    print('PASS: original-proof inequalities, actual-feedback history closure and full-trajectory local domain; not maximum ranges.')
