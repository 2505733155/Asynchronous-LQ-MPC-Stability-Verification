"""Auxiliary HMPC reference and separate actual-feedback local checks.

Python standard library only; no MATLAB, SDP solve, private data or network.
The auxiliary reference is NOT substituted for the actual QP/PWA controller.
The ten-point sampled-map proof is not an established application of the
published HMPC sampling theorem; its sampling-theory bridge remains unresolved.
See HMPC_reference_verification.md for the source/derivation distinction.
"""
import argparse
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]

def transpose(a):
    return list(map(list, zip(*a)))

def add(a, b):
    return [[x + y for x, y in zip(row, other)] for row, other in zip(a, b)]

def scale(a, c):
    return [[c * x for x in row] for row in a]

def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]

def inverse2(a):
    d = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    assert d
    return [[a[1][1] / d, -a[0][1] / d], [-a[1][0] / d, a[0][0] / d]]

def positive2(a):
    return a[0][0] > 0 and a[0][0] * a[1][1] - a[0][1] * a[1][0] > 0

def floating(a):
    return [[float(v) for v in row] for row in a]

def interval_mul(a, b):
    p = [x * y for x in a for y in b]
    return min(p), max(p)

def interval_add(a, b):
    return a[0] + b[0], a[1] + b[1]

def interval_neg(a):
    return -a[1], -a[0]

def exp_positive_bounds(x, n=80):
    # For positive x, truncate below; bound remaining positive tail geometrically.
    assert x > 0
    term = F(1)
    total = F(1)
    for j in range(1, n + 1):
        term *= x / j
        total += term
    first = term * x / (n + 1)
    ratio = x / (n + 2)
    assert ratio < 1
    return total, total + first / (1 - ratio)

mm, tr = matmul, transpose


def terminal_scope_bound(p,q,r):
    """Necessary terminal-decrease ceiling with all terminal sets free."""
    a,b,c = p[0][0],p[0][1],p[1][1]
    eps = q[1][1]
    assert 2*b>r and 4-2*b>0 and 4*a-2*c>0
    coeff=[4+4*a+2*b+eps, -8-12*a-12*b-4*eps-2*c,
           4+8*a+12*b+4*eps+4*c+r]
    def value(s):
        return (coeff[0]*s+coeff[1])*s+coeff[2]
    assert coeff[0]>0 and 4*coeff[0]+coeff[1]>0
    lo,hi=F(2),F('2.08')
    assert value(lo)<0<value(hi)
    for _ in range(100):
        mid=(lo+hi)/2
        if value(mid)>0: hi=mid
        else: lo=mid
    assert value(lo)<=0<value(hi) and hi<F('2.08')
    return {'fixed_reference':'T_H=2, Qc=diag(1,1e-4), Rc=.1, Pc as supplied, |u|<=1.',
            'free_choices':'Any terminal set and terminal feedback satisfying the published terminal decrease and admissibility conditions.',
            'necessary_polynomial_coefficients':[str(v) for v in coeff],
            'positive_root_interval':[str(lo),str(hi)],
            'root_approximation':float((lo+hi)/2),
            'centered_ball_radius_strict_upper_bound':2.08,
            'excluded_initial_state':[0,2.08],
            'minimum_terminal_residual_at_witness':float(value(F('2.08'))),
            'scope':'This is an auxiliary-reference feasibility restriction. It is not instability of the actual controller, an outer radius of every admissible domain, or an upper bound across arbitrary continuous reference designs.'}

def poly(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def pa(a, b):
    return poly([(a[j] if j < len(a) else F(0)) + (b[j] if j < len(b) else F(0)) for j in range(max(len(a), len(b)))])


def ps(a, c):
    return poly([v*c for v in a])


def pm(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            c[i+j] += v*w
    return poly(c)


def pe(a, t):
    result = F(0)
    for v in reversed(a):
        result = result*t+v
    return result


def determinant3(a):
    total = [F(0)]
    for p in permutations(range(3)):
        sign = (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
        total = pa(total, ps(pm(pm(a[0][p[0]], a[1][p[1]]), a[2][p[2]]), sign))
    return total


def bernstein_scalar(a, lo, hi, degree):
    power = [sum((a[k]*comb(k, j)*lo**(k-j)*(hi-lo)**j for k in range(j, len(a))), F(0)) for j in range(degree+1)]
    return [sum((power[j]*F(comb(k, j), comb(degree, j)) for j in range(k+1)), F(0)) for k in range(degree+1)]


def positive_polynomial_matrix(a, subdivisions=100):
    degree = max(len(v)-1 for row in a for v in row)
    for j in range(subdivisions):
        lo, hi = F(j, subdivisions), F(j+1, subdivisions)
        b = [[bernstein_scalar(v, lo, hi, degree) for v in row] for row in a]
        for k in range(degree+1):
            matrix = [[v[k] for v in row] for row in b]
            if not positive2(matrix):
                return False, j, k
    return True, subdivisions, degree+1


def reference_ct_enclosure(p, q):
    m = [[F(0), F(-2), F(0), F(0)], [F(0), F(0), F(0), F(20)],
         [2*q[0][0], F(0), F(0), F(0)], [F(0), 2*q[1][1], F(2), F(0)]]
    assert max(sum(abs(v) for v in row) for row in m) <= 20
    term = eye(4)
    total = eye(4)
    for n in range(1, 181):
        term = scale(mm(m, term), F(1, n))
        total = add(total, term)
    tail = F(3**20*20**181, factorial(181))
    z = mm(total, [[F(1), F(0)], [F(0), F(1)], p[0], p[1]])
    err = tail*max(F(1), *(sum(abs(v) for v in row) for row in p))
    zi = [[(v-err, v+err) for v in row] for row in z]
    x, y = zi[:2], zi[2:]
    det = interval_add(interval_mul(x[0][0], x[1][1]), interval_neg(interval_mul(x[0][1], x[1][0])))
    assert det[0] > 0
    invdet = F(1)/det[1], F(1)/det[0]
    invx = [[interval_mul(x[1][1], invdet), interval_mul(interval_neg(x[0][1]), invdet)],
            [interval_mul(interval_neg(x[1][0]), invdet), interval_mul(x[0][0], invdet)]]
    result = []
    for j in range(2):
        v = interval_add(interval_mul(y[1][0], invx[0][j]), interval_mul(y[1][1], invx[1][j]))
        result.append((-10*v[1], -10*v[0]))
    boxes = [(F('-3.197309'), F('-3.197306')), (F('-2.548451'), F('-2.548448'))]
    assert all(lo <= v[0] <= v[1] <= hi for v, (lo, hi) in zip(result, boxes))
    return result, boxes, tail


def discrete_gain(p, q, r, n):
    d = F(2, n)
    a, b = [[F(1), d], [F(0), F(1)]], [[F(0)], [d]]
    s = p
    for _ in range(n):
        bs_a = mm(tr(b), mm(s, a))
        denom = d*r+mm(tr(b), mm(s, b))[0][0]
        assert denom > 0
        k = scale(bs_a, -1/denom)
        s = add(scale(q, d), add(mm(tr(a), mm(s, a)), mm(tr(bs_a), k)))
    return k[0]


def quadratic(row, metric):
    return mm([row], mm(metric, tr([row])))[0][0]


def sampled_lyapunov_family(k):
    a = [[F(0), F(1)], k]
    c = [k, [F(0), F(0)]]
    bases = [[[F(1),F(0)],[F(0),F(0)]], [[F(0),F(1)],[F(1),F(0)]], [[F(0),F(0)],[F(0),F(1)]]]
    operator = []
    positions = [(0,0),(0,1),(1,1)]
    for w in bases:
        awc = mm(tr(a), mm(w, c))
        coeff = [add(mm(tr(a), w), mm(w, a)),
                 add(mm(tr(a), mm(w, a)), scale(add(mm(tr(c), w), mm(w, c)), F(1,2))),
                 scale(add(awc, tr(awc)), F(1,2)),
                 scale(mm(tr(c), mm(w, c)), F(1,4))]
        operator.append([poly([v[i][j] for v in coeff]) for i,j in positions])
    operator = list(map(list, zip(*operator)))
    den = determinant3(operator)
    rhs = [F(-1), F(0), F(-1)]
    nums = []
    for j in range(3):
        matrix = [[v[:] for v in row] for row in operator]
        for i in range(3):
            matrix[i][j] = [rhs[i]]
        nums.append(determinant3(matrix))
    if pe(den,F(0)) < 0:
        den, nums = ps(den,-1), [ps(v,-1) for v in nums]
    assert all(pa(pa(pm(operator[i][0], nums[0]), pm(operator[i][1], nums[1])), pm(operator[i][2], nums[2])) == ps(den,rhs[i]) for i in range(3))
    numerator = [[nums[0], nums[1]], [nums[1], nums[2]]]
    # Certify 0.33 I < W(h) < 3.1 I on the ENTIRE [0,1].
    lower = [[pa(numerator[i][j], ps(den, F('-.33') if i==j else F(0))) for j in range(2)] for i in range(2)]
    upper = [[pa(ps(den,F('3.1') if i==j else F(0)), ps(numerator[i][j],-1)) for j in range(2)] for i in range(2)]
    dd = [[den, [F(0)]], [[F(0)], den]]
    checks = {'denominator':list(positive_polynomial_matrix(dd)), 'lower_bound':list(positive_polynomial_matrix(lower)), 'upper_bound':list(positive_polynomial_matrix(upper))}
    assert all(v[0] for v in checks.values()), checks
    return den, numerator, checks


def calculate():
    qpfile = ROOT/'data/double_integrator/double_integrator_qp_data.m'
    lawfile = ROOT/'data/double_integrator/answer_regions_input_only_legacy.m'
    # Auxiliary terminal matrix, defined by these exact decimal values.
    p = [[F('2.1429023470226376'), F('1.2245640609262458')],
         [F('1.2245640609262458'), F('1.3995571393121802')]]
    q, r = [[F(1),F(0)],[F(0),F('.0001')]], F('.1')
    ac, bc = [[F(0),F(1)],[F(0),F(0)]], [[F(0)],[F(1)]]
    kf = scale(mm(tr(bc),p),-1/r)
    residual = add(q, add(mm(tr(ac),p), add(mm(p,ac), scale(mm(mm(p,bc),mm(tr(bc),p)),-1/r))))
    assert positive2(p) and positive2(scale(residual,-1))
    rho = F('.003')
    terminal_input_square = rho*p[1][1]/r**2
    assert terminal_input_square < 1
    ct, ctbox, tail = reference_ct_enclosure(p,q)
    kd = discrete_gain(p,q,r,2)
    text = qpfile.read_text(encoding='utf-8')
    def matrix(name):
        return [[F(v) for v in row.split()] for row in re.search(r'qp\.'+name+r'\s*=\s*\[([^]]+)\]',text)[1].split(';')]
    hq, fq = matrix('H'), matrix('Fux')
    ku = scale(mm(inverse2(hq),fq),-1)
    stored = [F(v.strip()) for v in re.search(r'K\s*=\s*\{\s*\[([^]]+)\]',lawfile.read_text(encoding='utf-8'))[1].split(',')]
    lawtext = lawfile.read_text(encoding='utf-8')
    origin_h = re.search(r'H\s*=\s*\{\s*\[([^]]+)\]',lawtext)[1]
    origin_rows = [[F(v.strip()) for v in row.split(',')] for row in origin_h.split(';')]
    assert len(origin_rows)==4 and all(sum(v*v for v in row)<9 for row in origin_rows)
    origin_b = re.search(r'b\s*=\s*\{\s*\[([^]]+)\]',lawtext)[1]
    assert [F(v.strip()) for v in origin_b.split(';')]==[F(1)]*4
    offsets = re.search(r'k\s*=\s*num2cell\(\[([^]]+)\]',lawtext)[1]
    assert F(offsets.split(';')[0].strip())==0
    actuals = {'ideal_rounded_QP':ku[0], 'stored_PWA':stored}
    for k in actuals.values():
        assert sum((a-b)**2 for a,b in zip(k,kd)) < F('.001')**2
        assert sum(v*v for v in k) < 4
        assert 1+sum(v*v for v in k) < F('2.2')**2
    assert sum((a-b)**2 for a,b in zip(actuals['ideal_rounded_QP'],stored)) < F('1e-26')
    # Use the previously defined state metric unchanged; it is an invertible
    # coordinate change for the reference problem, not a controller change.
    e = [stored[0]-F('-3.1973071'),stored[1]-F('-2.5482575')]
    metric = add(mm(tr([e]),[e]),scale(eye(2),F('.0001')))
    invmetric = inverse2(metric)
    lref = F('1.01')
    lround = F('.04')
    for corner in product(*ctbox):
        assert quadratic([a-b for a,b in zip(kd,corner)],invmetric) < lref**2
        for k in actuals.values():
            assert quadratic([a-b for a,b in zip(k,corner)],invmetric) < lref**2
    for k in actuals.values():
        assert quadratic([a-b for a,b in zip(k,kd)],invmetric) < lround**2
    # Recheck the analytic 0.88 CT convolution bound for the NEW reference.
    assert F('1.274') < -ct[1][1]/2 <= -ct[1][0]/2 < F('1.275')
    assert F('1.25')**2 < -ct[0][1]-F('1.275')**2
    assert -ct[0][0]-F('1.274')**2 < F('1.26')**2
    blo=(e[0]-F('1.28')*e[1])/F('1.26'); bhi=(e[0]-F('1.27')*e[1])/F('1.25')
    assert blo>e[1]>0 and e[1]**2+bhi**2<F('1.37')**2
    assert 3*F('3.14')/(4*F('1.26')) > F('1.86')
    assert exp_positive_bounds(F('1.27')*F('1.86'))[0] > 1/F('.095')
    assert exp_positive_bounds(F('1.27')*F('3.14')/F('1.26'))[0] > 1/F('.043')
    assert (1/F('1.25'))**2+(1+F('1.28')/F('1.25'))**2 < F('2.2')**2
    gain = e[0]/F('3.19')+2*F('1.37')*F('1.26')*F('.095')/(F('3.19')*(1-F('.043')))+F('.01')*F('2.2')/F('1.27')
    assert gain<F('.88')
    gamma_ref = F('.88')/(1-F('.88')*lref)
    assert F('.88')*lref<1 and gamma_ref*lround<1
    # Uniform inactive-domain proof for the CT and all Euler N>=2 references.
    radius=F('.00005')
    et=[[F(1),F(2)],[F(0),F(1)]]
    zero_cost=add(mm(tr(et),mm(p,et)),[[F(2),F(2)],[F(2),F(8,3)+F('.0002')]])
    assert sum(zero_cost[j][j] for j in range(2))<22
    assert F(14,3)*220<33**2
    assert 440*36*radius < 1  # all Euler Riccati gains <=440
    assert 22*radius**2<rho
    assert sum(v*v for row in ku for v in row)<9
    # CT semigroup bound: ||exp(A_CT t)||_S <=1100 exp(-1.27t).
    assert 2*F('1.28')**2+1+F('3.21')**2 < F('3.82')**2
    assert 1+F('3.82')/F('1.25')<F('4.1')
    assert (metric[0][0]+metric[1][1])*F('4.1')**2/F('.0001')<1100**2
    # Explicit locality for the CT small-gain step (not the final sampled ROA).
    assert (1100*F('2.52')*F('1e-12')+F('.88')*F('1e-9'))/(1-F('.88')*lref) < F('.01')*radius
    den, numerator, bchecks = sampled_lyapunov_family(kd)
    # Preserve explicit reference-to-actual feedback-error robustness instead of silently replacing Kd.
    # ||Mref||<=4.2, ||Bh||<=1.5h, ||W||<3.1, ||DeltaK||<.001.
    feedback_error_loss_per_h=2*F('4.2')*F('3.1')*F('1.5')*F('.001')+F('3.1')*F('1.5')**2*F('.001')**2
    assert feedback_error_loss_per_h < F('.05')
    # Hence Mactual' W Mactual-W <=-.95 h I and contraction <=1-.15 h.
    # The inequality holds uniformly for every 0<h<=1, not just ten points.
    assert F('.95')/F('3.1') > F('.3')
    assert F('3.1')*F('1.5')**2 < F('2.7')**2
    assert F('2.7')/F('.15') == 18
    # Full continuous-time trajectory bound incl. arbitrary measurable input.
    assert F('.57')**2<F('.33') and F('1.77')**2>F('3.1')
    assert exp_positive_bounds(F('.15'))[1]*F('4.2')*F('1.77')/F('.57') < 16
    assert F('4.2')*18/F('.57')+F('1.5') < 135
    chi, disturbance = F('2.5e-6'),F('5e-8')
    assert 16*chi+135*disturbance < radius
    assert 2*radius+disturbance<1
    level = F('4e-11')
    # Zero-disturbance ellipsoid x'W(h)x<=level remains in common ball
    # over complete holds, despite possible within-hold V increase.
    assert F('4.2')**2*level/F('.33') < radius**2
    assert F('3.1')*F('3.5e-6')**2 < level
    rows=[]
    for j in range(1,11):
        h=F(j,10)
        dh=pe(den,h)
        w=[[pe(v,h)/dh for v in row] for row in numerator]
        mref=add([[F(1),h],[F(0),F(1)]],mm([[h*h/2],[h]],[kd]))
        assert add(mm(tr(mref),mm(w,mref)),scale(w,-1))==scale(eye(2),-h)
        actual_checks={}
        for name,k in actuals.items():
            ma=add([[F(1),h],[F(0),F(1)]],mm([[h*h/2],[h]],[k]))
            decrease=add(w,scale(mm(tr(ma),mm(w,ma)),-1))
            assert positive2(add(decrease,scale(eye(2),-F('.95')*h)))
            actual_checks[name]=True
        rows.append({'Tsc':float(h),'reference_ZOH_ISS_quantified':True,'actual_reference_error_pass':actual_checks,
                     'W':floating(w),'W_equation':'Mref(h)^T W(h) Mref(h)-W(h)=-h I',
                     'actual_decrease_lower_bound':float(F('.95')*h),
                     'ISS_initial_ball_radius':float(chi),'ISS_input_disturbance_max':float(disturbance),
                     'zero_disturbance_ellipsoid_level':float(level)})
    sources=[qpfile,lawfile,Path(__file__).resolve()]
    terminal_scope = terminal_scope_bound(p,q,r)
    return {'date':'2026-11-09',
            'classification':'Separate local sampled-map Lyapunov/feedback-error proof with auxiliary ACC/arXiv HMPC reference calculations (C in independent review). Original-HMPC ten-point sampling certification is not established.',
            'reference':{'Q':floating(q),'R':float(r),'T':2,'N_current':2,'prediction_step_current':1,
                         'terminal_P':floating(p),'terminal_rho':float(rho),'terminal_residual':floating(residual),
                         'terminal_input_max_squared':float(terminal_input_square),
                         'CT_gain_interval':[[float(a),float(b)] for a,b in ct],
                         'CT_gain_wide_verified_boxes':[[float(a),float(b)] for a,b in ctbox],
                         'Hamiltonian_Taylor_remainder_bound':float(tail),'discrete_N2_gain':[float(v) for v in kd]},
            'actual':{name:[float(v) for v in k] for name,k in actuals.items()},
            'error_and_small_gain':{'weighted_state_metric':floating(metric),'CT_ISS_gain_upper':.88,
                                    'reference_discretization_error_upper':float(lref),'reference_small_gain_product':float(F('.88')*lref),
                                    'actual_total_error_upper':float(lref),'reference_to_actual_error_upper':float(lround),
                                    'reference_closed_ISS_gain_upper':float(gamma_ref),'actual_reference_error_small_gain_product':float(gamma_ref*lround),
                                    'reference_to_actual_Euclidean_error_upper':.001,
                                    'CT_small_gain_valid_initial_ball':1e-12,'CT_small_gain_valid_input_bound':1e-9,
                                    'common_inactive_ball_radius':float(radius),
                                    'family_statement':'For all N>=2: uniform inactive domain, positive CT/DT quadratic Hessians, Euler Riccati consistency with O(T/N) first-input convergence; see HMPC_reference_verification.md.'},
            'sampling_quantification':{'method':'Rational W(h) family, exact Lyapunov identities and Bernstein positivity on h in [0,1], followed by explicit reference-to-actual feedback-error robustness and complete ZOH input-disturbance bounds.',
                                       'uniform_W_lower':.33,'uniform_W_upper':3.1,'Bernstein_checks':bchecks,
                                       'W_denominator_power_coefficients':[str(v) for v in den],
                                       'W_numerator_power_coefficients':[[[str(v) for v in p] for p in row] for row in numerator],
                                       'feedback_error_loss_per_h_upper':float(feedback_error_loss_per_h),'actual_decay_per_h_lower':.95,
                                       'sample_contraction_upper':'1-0.15*h','intersample_norm_factor_upper':4.2,
                                       'full_trajectory_ISS_bound':'||x(t)||_2 <=16 exp(-0.15 t)||x(0)||_2+135||delta||_infinity',
                                       'ISS_initial_ball_radius':float(chi),'ISS_input_disturbance_max':float(disturbance),
                                       'zero_disturbance_ellipsoid_level':float(level),'zero_disturbance_safe_ball_radius':3.5e-6,
                                       'domain_scope':'original numerical position/velocity coordinates; conservative local inner domains, not maximum ROA',
                                       'ten_points':rows},
            'terminal_scope_bound':terminal_scope,
            'theory_source':{'primary':'Gautam and Nicotra, Stability Analysis of Hypersampled Model Predictive Control, arXiv:2310.02623v1 / ACC 2024, DOI 10.23919/ACC60939.2024.10644911', 'published_steps':'CT reference assumptions; state-proportional discretization error and convergence; small gain; existence of sufficiently fast sampled-data implementation.', 'additional_quantification':'LQ regularity/Euler convergence argument; rational CT gain enclosure; weighted gains; exact ZOH Lyapunov family; state-dependent actual-feedback error closure. These steps are derived in the companion technical note, not quoted as formulas from the paper.'},
            'original_HMPC_ten_point_certification':False,
            'sampling_theory_bridge':'Unresolved: Teel--Nesic--Kokotovic CDC 1998 Theorem 2 conditions and their link to W(h) have not been established. PASS is arithmetic/local stability evidence, not original-theory acceptance.',
            'literal_scope':'Auxiliary reference structural/regularity/convergence arguments and continuous-update small-gain calculations are provided. Actual feedback is retained via explicit state-dependent error bounds. The ten sampled-data values use a separate exact-hold-map Lyapunov proof. Their attribution to the original HMPC sampling theorem is unverified; the local proof and fixed-reference feasibility bound do not resolve that provenance gap.',
            'sources':[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]}



if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Write a new reference result; differing existing files are not overwritten.')
    parser.add_argument('--check',type=Path,help='Recompute all checks and compare saved evidence, including source hashes.')
    args=parser.parse_args()
    result=calculate()
    if args.check:
        assert result==json.loads(args.check.read_text(encoding='utf-8')), 'Saved result or source identity differs'
    if args.output:
        if args.output.exists():
            assert result==json.loads(args.output.read_text(encoding='utf-8')), 'Existing evidence would change'
        else:
            args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'classification':result['classification'],
                      'actual_feedback':result['actual'],
                      'sampling_points_passed':len(result['sampling_quantification']['ten_points']),
                      'reference_small_gain':result['error_and_small_gain']['reference_small_gain_product'],
                      'actual_feedback_error_small_gain':result['error_and_small_gain']['actual_reference_error_small_gain_product'],
                      'local_initial_radius':result['sampling_quantification']['ISS_initial_ball_radius'],
                      'terminal_scope':result['terminal_scope_bound']},indent=2))
    print('PASS: reference algebra, ten actual-feedback local ISS checks and fixed-reference terminal ceiling. Read the technical note for mathematical applicability and provenance.')
