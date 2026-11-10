# HMPC original-proof comparison and separate sampled-map verification

This note provides the quantitative HMPC comparison for the package's actual
QP/PWA controller. The original CT-reference, Euler-approximation and Teel
sampling proof, with explicit actual-feedback history closure, gives a local
certificate for every fixed 0<h<=.13 s and every initial state in the original
Euclidean closed ball of radius .0195, without external disturbance. The
fixed-reference derivation and both actual gain representations have been
independently checked. These are sufficient inner ranges, not maximum HMPC
capabilities. A separate sampled-map proof covers the ten tested values; its
ten-point result retains its distinct provenance.
The controller being verified is the package's existing rounded QP/PWA law;
the continuous and discrete reference controllers are auxiliary proof objects.

## Target and reproduction

From the package root, with Python 3 and its standard library:

```text
python -B src/hmpc/verify_hmpc_original_conditions.py --profile higher_period --check results/reference/hmpc_original_conditions_higher_period.json
python -B src/hmpc/verify_hmpc_reference.py --check results/reference/hmpc_reference_scope_results.json
```

The first entry checks the original-proof certificate; the second retains the
separate ten-point sampled-map calculation and terminal ceiling. Both use the
same controller data. The plant is `x_dot=[x2,u]`; the input is updated at a
fixed h and held constant until the next update. The paper grid is
`h in {0.1,0.2,...,1.0}`. Its exact map is

\[
x_{k+1}=\begin{bmatrix}1&h\\0&1\end{bmatrix}x_k+
\begin{bmatrix}h^2/2\\h\end{bmatrix}\kappa(x_k).
\]

The actual controller reads `data/double_integrator/double_integrator_qp_data.m`:

\[
\min_{|U_j|\le1}\;\tfrac12U^\top HU+U^\top F_{ux}x,\quad
H=\begin{bmatrix}.8365&.3603\\.3603&.2059\end{bmatrix},\quad
F_{ux}=\begin{bmatrix}.4624&1.2852\\.1682&.5285\end{bmatrix}.
\]

There is no terminal-set constraint in this actual QP. The archived table in
`answer_regions_input_only_legacy.m` has origin feedback
`K0=[-.815801595784556,-1.749302793118031]`. The ideal rounded-QP origin gain
is the first row of `-H^-1 Fux`. Their difference is checked below `1e-13`;
they are verified separately rather than treated as identical finite decimals.
The first region's inequalities and zero affine offset are also read from the
table. No gain from an auxiliary controller replaces either actual gain.

All radii below use the original numerical position/velocity coordinates.
They are not position-only distances or radii in the manuscript's SDP scaling.

## Published theory and the quantitative steps supplied here

Primary source: Y. Gautam and M. M. Nicotra, *Stability Analysis of Hypersampled
Model Predictive Control*, [arXiv:2310.02623v1](https://arxiv.org/abs/2310.02623),
published at ACC 2024, DOI `10.23919/ACC60939.2024.10644911`.
The ACC terminal inequality includes the derivative along the dynamics,
`grad(J) f+l<=0`; the local preprint's missing dynamics factor is a typesetting
difference, not a different terminal condition.

| Source location | Published role | Instance or additional derivation here |
| --- | --- | --- |
| Assumption 1, Theorem 1 | Stable CT-OCP with terminal conditions; local ISS | Positive quadratic reference cost, feasible terminal ellipsoid, local LQ solution and a computed input/state gain |
| Assumption 2 | CT/DT strong regularity and state-proportional feedback error tending to zero under refinement | Constraint-inactive LQ neighborhood, genuine Euler family, Riccati consistency argument below |
| Theorem 2 and its small-gain inequality | Replace CT feedback by the discrete-OCP feedback | Same weighted norm for error and ISS gain; fixed-step product `0.8888<1` |
| Theorem 3 | Sufficiently fast ZOH implementation preserves local ISS, via Teel--Nesic--Kokotovic (CDC 1998), Theorem 2 | The original inequalities are quantified below; after actual-feedback closure they cover fixed h<=.13, hence .1 on the paper grid. The W(h) ten-point proof remains separate |
| Input-disturbance interface | Admits an additive input disturbance | The original-proof entry closes the actual error in the same input position, with a 3h state history and the original Teel Theorem 1; both actual gains are checked |

The reference problem and its Euler family satisfy the structural/regularity
conditions used by ACC Theorems 1--2. The sampling proof below explicitly
instantiates Teel--Nesic--Kokotovic, CDC 1998, Theorems 1--2, (10)--(18) and
Comment 1 (DOI `10.1109/CDC.1998.757793`). The original sampling conditions
are checked, rather than inferred from a stable exact hold map. The actual
QP is not asserted identical to the reference OCP: its state-dependent input
error is closed by the original history criterion. The master's-thesis
alternative ordering is not used. The independent proof review accepted this
as a quantitative instantiation with explicit feedback-error treatment.

## Auxiliary reference and a common local domain

The auxiliary CT-OCP has physical horizon `T_H=2`, input bound `|u|<=1`, and

\[
Q_c=\operatorname{diag}(1,10^{-4}),\quad R_c=.1,\quad
P_c=\begin{bmatrix}
2.1429023470226376&1.2245640609262458\\
1.2245640609262458&1.3995571393121802
\end{bmatrix}.
\]

These displayed decimals define the reference matrix. Its terminal set is
`x^T Pc x<=.003`. With `Kf=-Rc^-1 Bc^T Pc`, rational Sylvester tests give
`Pc>0` and

\[
A_c^\top P_c+P_cA_c-P_cB_cR_c^{-1}B_c^\top P_c+Q_c<0.
\]

The squared maximum terminal input is below `.419868<1`; the terminal
ellipsoid is forward invariant. The small positive velocity weight ensures
the reference cost has a strict state-quadratic lower bound. It does not
change the actual QP's weights, matrices or constraints.

For every `N>=2`, define the actual Euler approximation of this one CT-OCP by
`d=2/N`, `Ad=I+d Ac`, `Bd=d Bc`, `Qd=d Qc`, `Rd=d Rc`, with the same `Pc`
and terminal set. This is the reference family, not a claim that the actual
rounded QP has those unrounded condensed coefficients.

The ball `||x||_2<=r=5e-5` is a uniform constraint-inactive domain:

* The zero-input prediction cost matrix has trace below 22 for CT and all
  these Euler grids. For the double integrator `Ad^j=I+jd Ac`; left sums of
  the positive increasing diagonal terms are bounded by their integrals.
* The unconstrained optimum therefore has cost at most `22||x||^2` and
  control energy at most `220||x||^2`. Variation of constants and
  Cauchy--Schwarz give predicted state norm below `36||x||`.
* The Riccati matrices have norm at most 22. CT gains have norm at most 220;
  Euler gains satisfy `||(Rc+d Bc^T S Bc)^-1 Bc^T S(I+d Ac)||<=440`.
  All predicted inputs are at most `440*36*r=.792<1`; terminal cost is
  below `22*r^2<.003`.

The actual QP's two unconstrained input rows have combined norm below 3.
The stored origin-region rows have norm below 3 and right sides one.
Both actual descriptions are linear on the same ball.

For local strong regularity, eliminate the linear state equation. The CT
reduced input Hessian is `2 Rc I` plus bounded positive semidefinite integral
and terminal operators, hence coercive and invertible on `L2`. For bounded
data perturbations its inverse solution is also bounded: the optimality
equation is a positive multiple of the input plus an integral term mapping
`L2` inputs into bounded functions. The strict constraint slack then keeps
the inequalities inactive in a sufficiently small `L-infinity` neighborhood.
The linear state equation is solvable with a bounded inverse. The resulting
local KKT system is invertible and the solution is single-valued and
Lipschitz. The DT reduced Hessian is at least `2d Rc I>0`, with full-row-rank
state equalities and the same inactivity argument. This analytic argument,
not just a test of the actual two-variable Hessian, addresses the CT/DT
regularity assumption.

For feedback convergence, the Euler Riccati map is

\[
\Phi_d(S)=dQ_c+(I+dA_c)^\top
\left[S-dSB_c(R_c+dB_c^\top SB_c)^{-1}B_c^\top S\right](I+dA_c).
\]

On the common bounded positive-semidefinite matrix set it equals
`S+d F(S)+O(d^2)` uniformly, where
`F(S)=Ac^T S+S Ac-S Bc Rc^-1 Bc^T S+Qc` is the CT Riccati field.
Its derivative is `(Ad+Bd K(S))^T E(Ad+Bd K(S))`; the bound on `K(S)` gives
a one-step Lipschitz constant `1+O(d)`. Discrete Gronwall on the fixed horizon
therefore gives Riccati and first-input errors `O(d)`. This proves existence
of an error function tending to zero on the common neighborhood, rather than
inferring a limit from a finite refinement plot. No monotonic error envelope
from the master's thesis is assumed.

## Fixed prediction step and feedback errors

For `d=1,N=2`, the reference DT gain is approximately
`Kd=[-.816562685446818,-1.7498807770434979]`.
The CT finite-horizon gain is enclosed by
`Kc1 in [-3.197309,-3.197306]`, `Kc2 in [-2.548451,-2.548448]`.
The latter follows from a 180-term rational Hamiltonian exponential with a
proved remainder bound. Float displays of tiny interval widths are not
used to decide signs.

Let `e=K0-[-3.1973071,-2.5482575]` and `S=e^T e+1e-4 I`.
In the norm `||x||S=sqrt(x^T Sx)`, interval corner checks give
`|(Kd-Kc)x|<=1.01||x||S`, and each actual/reference error satisfies
`|(Kactual-Kd)x|<=.04||x||S` and `||Kactual-Kd||_2<.001`.

The CT impulse response is

\[
g(t)=e^{-at}[\sin(bt)/b,\;\cos(bt)-(a/b)\sin(bt)]^\top,
\quad a=-K_{c2}/2,\quad b=\sqrt{-K_{c1}-a^2}.
\]

The checker bounds `1.27<a<1.28`, `1.25<b<1.26`. Using
`||g||S<=|eg|+.01||g||2`, the absolute integral of the scalar damped sinusoid,
and `||g(t)||2<=2.2 exp(-at)`, it obtains an integral bound below `.88`.
Thus the fixed-step small-gain product is `.88*1.01=.8888<1`.
The reference continuously updated DT-feedback asymptotic gain is bounded by
`.88/(1-.8888)=7.913669065`; multiplying by the actual feedback-error bound
`.04` gives `.316547<1`. The error vanishes at the origin.

Local validity is retained, not inferred from an asymptotic gain alone.
The CT semigroup has weighted bound `1100 exp(-1.27t)`; finite-time supremum
and small-gain bounds provide an initial neighborhood and disturbance limit
remaining in the legacy common ball. The W(h) branch below retains those
small local bounds. The original-proof branch next gives a larger common
domain directly from the fixed reference constraints.

## Quantitative original sampling proof and actual-feedback closure

The portable entry `verify_hmpc_original_conditions.py` uses standard-library
rational arithmetic. It reads the actual QP and stored PWA data, recomputes
the fixed reference quantities and checks all three profiles separately:

| Profile | Every fixed period h / s | Initial Euclidean ball radius | Input disturbance norm bound | Reference input gain Gamma |
| --- | --- | --- | --- | --- |
| higher_disturbance | 0<h<=.11 | .0195 | .0005 | 22 |
| balanced | 0<h<=.12 | same | <=.00025 | 48 |
| higher_period | 0<h<=.13 | same | <=.00003 | 377 |

Only h=.1 belongs to both these intervals and the paper's ten-point grid.
The period is fixed in each implementation; arbitrary switching between
periods is not asserted. Each statement covers all initial states in the
closed ball and the complete continuous trajectory, including between updates.
For w=0 the origin is locally asymptotically stable. The larger period profile
allows less disturbance, so the rows are not equal-robustness comparisons.

### Original Lyapunov and sampling conditions

For the reference Kd=(-a,-b), take

\[
Y=\begin{bmatrix}1&.98\\.98&1.4272\end{bmatrix},\qquad
V=x^\top Yx,\qquad s=\|x\|_Y.
\]

The change of coordinates z=Lx, L^T L=Y, makes the comparison functions
alpha1(s)=alpha2(s)=s^2. All gains below use this same norm. Set
DY=-[(Ac+BcKd)^T Y+Y(Ac+BcKd)], q<=lambda_min(DY,Y),
p>=sqrt(Bc^T Y Bc), and

\[
C\ge\|K_dA_c\|_{Y^{-1}}+b\|K_d\|_{Y^{-1}}.
\]

The rational witness has q>=1.26964092318649, p<=1.194654761846 and
C<=4.016060054380895. With theta=Kd[x(kh)-x(t)] and additive input delta3,

\[
\dot V\le-qs^2+2ps(|\theta|+|\delta_3|),\quad
|\theta|\le h\max\{(1+\varepsilon)CX_{2h},
(1+\varepsilon^{-1})bD_{2h}\}.
\]

X denotes the state-history Y norm and D the input-history essential supremum.
For positive allocations rho_theta+rho_input<1 choose

\[
c_\theta=\frac{2p}{q\rho_\theta},\quad
c_\delta=\frac{2p}{q\rho_{\mathrm{input}}},\quad
\alpha_3(s)=2q(1-\rho_\theta-\rho_{\mathrm{input}})s^2.
\]

The checker verifies both the continuous ISS premise and the original (15),
Vdot<=-alpha3/2 at the two input thresholds. With G=1, original (18) becomes

\[
g_h=c_\theta h(1+\varepsilon)C<1,\qquad
\Gamma_h\ge\max\{c_\theta h(1+\varepsilon^{-1})b,c_\delta\}.
\]

The three allocations (rho_theta,rho_input,epsilon) are
(.85,.10,.02), (.94,.04,.01), and (.99,.005,.002). At their period upper
endpoints gh is respectively below .997622, .974466, and .994415. These
positive linear envelopes also cover every smaller positive h. The
higher_period witness's strict sufficient threshold is about .13073015;
.13 is an inward choice, not a maximum-period computation. The global linear
extensions satisfy the strict contraction for every s>0, permitting the
zero-offset Comment 1 case; the domain argument below returns to actual laws.

### State-dependent actual input and the original 3h history condition

For each actual representation separately the checker verifies

\[
e(x)=(K_a-K_d)x,\qquad |e(x)|\le\ell\|x\|_Y,\quad\ell=.0008.
\]

Write delta3(t)=e(x(kh))+w(t). Its original 2h input history refers to sampled
states as early as t-3h. Therefore

\[
D_{2h}\le\ell X_{3h}+W_{3h}
\le\max\{2\ell X_{3h},2W_{3h}\}.
\]

Substitution into original (17) gives

\[
s(t)\ge\max\{r_hX_{3h},2\Gamma_hW_{3h}\}
\ \Longrightarrow\ \dot V\le-\alpha_3(s(t))/2,\qquad
r_h=\max\{g_h,2\Gamma_h\ell\}<1.
\]

This directly matches the original Teel Theorem 1, whose history length is
arbitrary: here it is 3h. Its functions are gamma_x(s)=r_h s and
gamma_w(s)=2 Gamma_h s, with alpha1=alpha2=s^2. The additional feedback gains
2 Gamma_h ell are .0352, .0768 and .6032. This closes the state-dependent error;
merely treating it as a bounded external disturbance would not give convergence.
The remaining external signal is w. No W(h) equation or different sampling
theorem is used. Physical startup uses constant state prehistory and the first
integral starts at zero; bounded measurable w needs no time derivative.

### Fixed-reference valid domain and the complete trajectory

Set Pbar=1.01 Pc as a proof upper bound, not a changed terminal cost. With the
stored origin row K0 as an unconstrained bounding policy and A0=Ac+BcK0,

\[
\Delta_1=\bar P-(I+A_0)^\top\bar P(I+A_0)-Q_c-R_cK_0^\top K_0\succ0.
\]

For every 0<d<=1 the Euler residual is

\[
\bar P-(I+dA_0)^\top\bar P(I+dA_0)-d(Q_c+R_cK_0^\top K_0)
=d[\Delta_1+(1-d)A_0^\top\bar PA_0]\succ0.
\]

The continuous residual is negative as well. Telescoping/integration bounds
all unconstrained CT and Euler N>=2 value matrices by Pbar. Along optimal
predictions the remaining cost and terminal cost do not exceed x0^T Pbar x0.
The CT stationarity equation u=-Rc^-1 Bc^T Sx and the Euler version
u_j=-Rc^-1 Bc^T S_(j+1)x_(j+1) imply

\[
|u|^2\ \text{or}\ |u_j|^2
\le\bar P_{22}R_c^{-2}x_0^\top\bar Px_0.
\]

Thus x^T Pbar x<.003 ensures all predicted input squares are below
.424065814<1 and the terminal cost is strictly below .003, for the CT whole
horizon and every Euler N>=2. Unconstrained optima therefore satisfy the actual
reference constraints. This is a uniform family argument, not finite N sampling.
The reduced CT Hessian is coercive on L2; continuous control/stationarity spaces
and the strictly inactive branch supply the local invertible stationarity map.
No open interior of a pointwise input box in L2 is assumed. The bounded Riccati
consistency argument above applies on this larger common domain.

The checker proves Pbar<3.522Y and chooses E={x:x^T Yx<=.00085}.
Since 3.522*.00085=.0029937<.003, E lies strictly inside that valid domain.
Each actual QP input row and PWA origin-region face is checked on E. The
original strict history criterion yields the all-time maximum bound

\[
\|x(t)\|_Y\le\max\{\|x_0\|_Y,2\Gamma_h\|w\|_\infty\}.
\]

With lambda_max(Y)<2.216608, (.0195)^2 lambda_max(Y)<.00085. Each tabulated
external-input bound also has (2 Gamma_h wbar)^2<.00085. A first-exit argument
therefore makes the linear extension and actual local controller coincide for
all time; constraints and intersample trajectories are covered. Attraction
follows from the same original history theorem, and for w=0 the state tends to
zero. The bound does not assert that V decreases throughout each hold interval.

All these radii are sufficient inner estimates for the fixed reference and
proof choices. A fixed-terminal-set feasibility bound gives
.0195<=r*(h)<.720081 for 0<h<=.13; r* is the supremum of full centered balls
certifiable through this reference. Allowing terminal sets/feedbacks to change
while holding the stated cost, horizon and input bound fixed instead retains
the 2.08 necessary ceiling below. Neither upper bound is an outer actual ROA
or a universal limit over redesigned HMPC references. The quadratic/gain
choices have not been proved globally optimal.

## Separate sampled-map proof and retention of the actual feedback

Set `Ar=Ac+Bc Kd`, `Cr=Ac Bc Kd`, and
`Mr(h)=I+h Ar+(h^2/2)Cr`. The checker constructs a rational matrix function
`W(h)` satisfying

\[
M_r(h)^\top W(h)M_r(h)-W(h)=-hI,\qquad .33I<W(h)<3.1I.
\]

The equation is divided by `h` and extended to zero before solving the
three-entry symmetric linear system by rational Cramer formulas. Denominator
positivity and both matrix bounds are proved by Bernstein coefficients on
100 subintervals of `[0,1]`. This is a sufficient condition for stability of
the exact hold map. No implication from this condition to the hypotheses or
sampling threshold of the cited HMPC/Teel result has been established here.

For either actual gain write `DeltaK=Kactual-Kd`. Then
`Ma=Mr+Bh DeltaK`, `Bh=[h^2/2,h]^T`. The extra Lyapunov term is bounded by

\[
2\|M_r\|\|W\|\|B_h\|\|\Delta K\|
+\|W\|\|B_h\|^2\|\Delta K\|^2
\le .039066975h<.05h.
\]

Here `||Mr||<=4.2`, `||Bh||<=1.5h`, `||W||<3.1`, and `||DeltaK||<.001`.
Consequently `Ma^T W Ma-W<=-.95h I`. All ten points and both actual
representations are separately rechecked by exact rational identities.
This closes the state-dependent feedback error within the separate sampled-map
proof. It has the same additive-input form as the published disturbance
interface, but that fact alone does not attribute this W-based ten-point
argument to the original sampling theorem.
Treating the error as merely a bounded exogenous signal would not establish
convergence to the origin.

For an additional arbitrary bounded measurable input disturbance `delta(t)`,
the one-step integral increment has norm at most `1.5h||delta||infinity`.
With `zk=sqrt(xk^T W(h)xk)` this gives
`z[k+1]<=(1-.15h)zk+2.7h||delta||infinity`, hence
`zk<=exp(-.15kh)z0+18||delta||infinity`.
The hold-interval map has norm at most 4.2. Thus the complete trajectory obeys

\[
\|x(t)\|_2\le16e^{-.15t}\|x(0)\|_2+135\|\delta\|_\infty.
\]

For `||x0||2<=2.5e-6`, `||delta||infinity<=5e-8`, the bound stays below
`4.675e-5<r`, closing a first-exit argument. With zero disturbance the
ellipsoids `x^T W(h)x<=4e-11` also remain inside the common ball over complete
holds and contain radius `3.5e-6` balls. These conservative inner estimates
are not maximum attraction domains. Each `h` is fixed; arbitrary switching
between sampling values is not asserted.

## A terminal-feasibility ceiling independent of terminal-set tuning

Keep `T_H=2,Qc,Rc,Pc` and `|u|<=1`, but allow any terminal set and any
terminal feedback satisfying the original terminal conditions. Write
`Pc=[[a,b],[b,c]]`, `epsilon=1e-4`. Its required terminal residual is

\[
g(z,u)=z_1^2+2az_1z_2+(2b+\epsilon)z_2^2
+2(bz_1+cz_2)u+R_cu^2\le0.
\]

Starting from `x0=(0,s)`, `s>=2`, every reachable terminal state has
`z1>=2s-2>=2`, `z2>=s-2>=0`. On this quadrant the minimum residual over
`|u|<=1` occurs at `u=-1`, and it increases in both state coordinates.
Therefore every terminal residual is at least

\[
G(s)=(4+4a+2b+\epsilon)s^2
-(8+12a+12b+4\epsilon+2c)s
+4+8a+12b+4\epsilon+4c+R_c.
\]

Rational arithmetic gives `G(2.08)>0.007816` and its larger root below 2.08.
Thus no admissible terminal set/feedback for these fixed reference data can
cover `x0=(0,2.08)`. In particular, it cannot certify a complete centered
ball of radius 2.08. This statement neither excludes all states with norm
above 2.08 nor bounds every possible HMPC reference design. It is a
reference-feasibility restriction, not instability of the actual QP feedback.

## Evidence interpretation

The original-proof result files record all profiles, both actual gains,
strict sampling/history gains, the common reference-domain residual and the
original-coordinate initial ball, together with hashes of the package inputs.
The portable arithmetic reproduces the independently checked inequalities.
It preserves the actual/reference difference instead of claiming equality of
OCPs. The reported period/ball are sufficient ranges for the chosen reference,
not maxima; only h=.1 on the paper grid is covered by this proof witness.

The unchanged legacy result contains the separate W(h) family, ten checks and
terminal-residual polynomial. Its PASS remains local sampled-map arithmetic
and the terminal ceiling; it is not an original-theorem certificate for all ten
periods. The new original-proof entry does not upgrade that legacy argument.
Neither entry invokes a trajectory or SDP run.
