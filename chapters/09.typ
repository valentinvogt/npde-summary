#import "../src/setup.typ": *
#show: thmrules

= Second-Order Linear Evolution Problems
<ch:9>
#counter(heading).update((9,1))
== Parabolic Initial-Boundary Value Problems
<sub:parabolic-ivps>
=== Heat Equation

In local form the heat equation is given by
#neq($ frac(partial, partial t) (rho u) - div (kappa (bx) grad u) = f quad upright("in ") tilde(Omega) = Omega times openint(0, T) $) <eq:heat-local>
where $u$ is the temperature, $rho$ the heat capacity, $kappa$ the heat
conductivity and $f$ a (time-dependent) heat source/sink. Without the
time derivative, this looks very similar to standard PDEs, for which we
know the transformation into a nice variational problem.

To solve it we still need boundary conditions. Besides the boundary
conditions of the spatial domain --- which are now required for all times --- one also needs initial conditions over the whole domain at time 0.

#neq($ u (bx , t) & = g (bx , t) quad upright("for ") (bx , t) in partial Omega times openint(0, T)\
u (bx , 0) & = u_0 (bx) quad upright("for all ") bx in Omega $) <eq:heat-bc-ic>

Testing with time-independent test functions $v$ and assuming $rho$ to
be time independent as well, we get to
#neq($ integral_Omega rho (bx) dot(u) v dif bx + integral_Omega kappa (bx) med grad u dot.op grad v dif bx = integral_Omega f (bx , t) thin v dif bx quad forall v in H_0^1 (Omega)\
u (bx , 0) = u_0 (bx) in H_0^1 (Omega) $) <eq:heat-integral-form>

with the shorthand notation
$ m (dot(u) , v) & = integral_Omega rho (bx) med dot(u) med v dif bx\
a (u , v) & = integral_Omega kappa (bx) med grad u dot.op grad v dif bx\
ell (v) & = integral_Omega f (bx , t) med v dif bx $
and the realisation that $m (dot(u) , v) = frac(dif, dif t) m (u , v)$
as only $u$ depends on time (note that importantly, the domain $Omega$ also stays constant) we can rewrite @eq:heat-integral-form as
#neq($ frac(dif, dif t) m (u , v) + a (u , v) = ell (v) $) <eq:heat-short>

which looks like something we know how to solve from NumCSE.

#pagebreak()
#counter(heading).update((9,2,2))
=== Stability

#lemma(number: "9.2.3.8.", "Decay of solutions of parabolic evolutions")[
  $colMath("If " f equiv 0, accentcolor)$, the solution $u (t)$ of the heat equation @eq:heat-integral-form satisfies
  $ norm(u (t))_m lt.eq e^(- gamma t) norm(u_0)_m , quad norm(u (t))_a lt.eq e^(- gamma t) norm(u_0)_a quad forall t in openint(0, T) $
  where $gamma = upright("diam") (Omega)^(- 2)$.
]
Note that this lemma also tells us that if $f$ is time-independent, the
solution $u (t)$ converges exponentially (in time) to the stationary
solution (the solution of @eq:heat-short without the
$m (dot.op , dot.op)$ part).

=== Method of Lines

Now let's look into how we can solve @eq:heat-short. Let's apply the Galerkin discretization. As $u$ is now
also time-dependent, we let the coefficients of $u$ (and not the basis) depend on time. $ u_h (t) = sum_(i = 1)^N mu_i (t) b_h^i $ Combining this
with @eq:heat-short, we get
#neq($ bold(M) {frac(dif, dif t) arrow(mu) (t)} + bA arrow(mu) (t) & = arrow(phi) (t)\
arrow(mu) (0) & = arrow(mu)_0 $) <eq:heat-galerkin>

where $bold(M)_(i , j) = m (b_h^j , b_h^i)$, $bA_(i , j) = a (b_h^j , b_h^i)$ and $[arrow(phi) (t)]_i = ell (b_h^i)$.

This is now an ODE with respect to time and can be solved by time
stepping, learned in NumCSE.

#strong[Recall ODEs]
An ODE is given as
$ bold(dot(u)) = bold(f \() t , bold(u) \) $ and is called linear if
$bold(f \() t , bold(u) \) = bA (t) bold(u)$. There is an evolution operator associated with the ODE, defined as
$Phi^(t_0 , t) u_0 = u (t)$. There are some methods to approximate the
evaluation operator with a discrete evolution operator $Psi$.

- explicit Euler:
  $Psi^(t , t + tau) bu = bu + tau bold(f) \( t , bu \)$

- implicit
  Euler:$Psi^(t , t + tau) bu = bw \, bw = bu + tau bold(f) \( t + tau , bold(w) \)$

- implicit midpoint:
  $Psi^(t , t + tau) bu = bw \, bw = bu + tau bold(f) \( t + 1 / 2 tau , 1 / 2 (bold(w + u \)))$

Hence we can calculate the time evolution by the sequence
$ bold(u)^((0)) = bold(u)_0 , quad bold(u)^((j)) = Psi^(t_(j - 1) , t_j) bold(u)^((j - 1)) , quad j = 1 , dots.h , M $

As $Psi$ is the discrete approximation, the question about the error is
immediate. One usually considers

- the error at final time:
  $epsilon.alt_M = norm(bold(u) (T) - bold(u)^((M)))$

- maximum error in the sequence:
  $epsilon.alt_oo = max_j norm(bold(u)^((j)) - bold(u) (t_j))$

#theorem(number: "9.2.6.14", "Convergence of single-step methods")[
  Given the above sequence of solutions, obtained by a single step method of order $q in bb(N)$, then
  $ epsilon.alt_oo = max_j norm(bold(u)^((j)) - bold(u) (t_j)) lt.eq C tau^q $
  with $tau = max_j lr(|t_j - t_(j - 1)|)$.
]

#strong[Runge--Kutta Single-Step Methods]
#definition(number: "7.3.3.1", [General Runge--Kutta single-step method])[
  For coefficients $b_i , a_(i , j) in bb(R) , c_i = sum_(j = 1)^s a_(i , j)$, the discrete evolution operator $Psi^(s , t)$ of an #strong[s-stage Runge--Kutta single step method] (RK-SSM) for the ODE $bold(dot(u)) = bold(f \() t , bold(u) \)$ is defined by
  $ bold(k)_i = bold(f) (t + c_i tau , bold(u) + tau sum_(j = 1)^s a_(i , j) bold(k)_j) , quad i = 1 , dots.h , s , quad Psi^(t , t + tau) bold(u) = bold(u) + tau sum_(j = 1)^s b_j bold(k)_j $
  with $bold(k)_j$ the increments.
]

The RK-SSM methods can be written down in compact form (the butcher
scheme) as 
#set math.mat(gap: 1em)
#neq($ mat(
  bold(c), bold(frak(A));  "", bold(b);
  delim: #none, 
  augment: #(
  hline: 1,
  vline: 1,
  stroke: 1pt
))
$) <eq:butcher>
 where $bold(c)$ is a vector containing the
coefficients $c_i$, $bold(b)$ the coefficients $b_i$ and $bold(frak(A))$ a
matrix containing the coefficients $a_(i , j)$.

So continuing from @eq:heat-galerkin with different time stepping
schemes, we get

- explicit Euler:
  $ arrow(mu)^((j)) = arrow(mu)^((j - 1)) + tau_j bold(M)^(- 1) (arrow(phi) (t_(j - 1)) - bA arrow(mu)^((j - 1))) $

- implicit Euler:
  $ arrow(mu)^((j)) = (tau_j bA + bold(M))^(- 1) (bold(M) arrow(mu)^((j - 1)) + tau_j arrow(phi) (t_(j - 1))) $

- implicit midpoint (Crank-Nicolson)
  $ arrow(mu)^((j)) = (bold(M) + 1 / 2 bA)^(- 1) tau_j ((bold(M) - 1 / 2 bA) arrow(mu)^((j - 1)) + 1 / 2 (arrow(phi) (t_j) + arrow(phi) (t_(j - 1)))) $

These all involve solving a linear system of equations each time step.
However note, that the matrices to invert stay constant with respect to
time, so we can calculate the decomposition only once to save a lot of
time.

Using a general RK-SSM method as the time step, we get the following
system of equations
$ bold(M) arrow(kappa)_i + sum_(m = 1)^s tau a_(i , m) bA arrow(kappa)_m & = arrow(phi) (t_j + c_i tau) - bA arrow(mu)^((j))\
arrow(mu)^((j + 1)) & = arrow(mu)^((j)) + tau sum_(m = 1)^s b_m arrow(kappa)_m $
With the Kronecker product, this can be rewritten as
$ (bold(I)_s times.circle bold(M) + tau bold(frak(A)) times.circle bA) mat(delim: "[", arrow(kappa)_1; dots.v; arrow(kappa)_s) = mat(delim: "[", arrow(phi) (t_j + c_1 tau) - bA arrow(mu)^((j)); dots.v; arrow(phi) (t_j + c_s tau) - bA arrow(mu)^((j))) $
which can be used to solve for the increments $arrow(kappa)_i$.

Recall stiff initial value problems:

#mybox("Stiffness")[
  An initial value problem is called stiff if stability imposes much tighter timestep constraints on explicit single step methods than the accuracy requirements.
]

To study the stiffness of the method of lines, we first diagonalize it.
For @eq:heat-galerkin let $arrow(psi)_1 , dots.h , arrow(psi)_N$ denote
the $N$ linearly independent generalized eigenvectors satisfying
$ bA arrow(psi)_i = lambda_i bold(M) arrow(psi)_i , quad (arrow(psi)_j)^top bold(M) arrow(psi)_i = delta_(i j) $
with positive eigenvalues $lambda_i$. With
$bold(T) = [arrow(psi)_1 , dots.h , arrow(psi)_N]$ and
$bold(D) = upright("diag") (lambda_1 , dots.h , lambda_N)$, this can be
rewritten as $ bold(A T = M T D \, quad T^top M T = I) $ The existence of eigenvectors with positive eigenvalues is guaranteed, as $bold(A \, M)$
are positive (semi)definite. Thus with a change of basis to the
eigenvector basis, one can diagonalize @eq:heat-galerkin.
$ arrow(mu) (t) = sum_k eta_k (t) arrow(psi)_k arrow.l.r.double arrow(mu) (t) = bold(T) arrow(eta) (t) arrow.l.r.double bold(T^top M) arrow(mu) (t) = arrow(eta) (t)\
arrow.r bold(M T) frac(dif, dif t) arrow(eta) (t) + bold(M T D) arrow(eta) (t) = arrow(phi) (t)\
arrow.r frac(dif, dif t) arrow(eta) (t) + bold(D) arrow(eta) (t) = bold(T)^top arrow(phi) (t) $
As $bold(D)$ is diagonal, this amounts to $N$ decoupled scalar ODEs. On those, we can perform our analysis more easily. In NumCSE you have seen that both Euler and Crank-Nicolson can be rewritten as a
RK-SSM with appropriate coefficients, so we can study the stability of
the general RK-SSM for the scalar case. For $dot(u) = - lambda u$, with
the butcher scheme @eq:butcher we obtain
$Psi_lambda^(t , t + tau) u = S (- lambda tau) u$, with the stability
function
$ S (z) = 1 + z bold(b)^top (I - z bold(frak(A)))^(- 1) bold(1) = frac(upright("det") (bold(I) - z bold(frak(A)) + z bold(b 1)^top), upright("det") (bold(I) - z bold(frak(A)))) $

#strong[Unconditional stability of single step methods] A necessary
condition for unconditional stability of a single step method, is that
the discrete evolution operator $Psi_lambda^t$ applied to the scalar ODE
$dot(u) = - lambda u$ satisfies $ lambda > 0 arrow.r lim_(j arrow.r oo) (Psi_lambda^tau)^j u_0 = 0 , quad forall u_0 , forall tau > 0 $
#definition(number: "9.2.7.46", "L-stability")[
  An RK-SSM satisfying the above condition, is called L-stable if its
  stability function satisfies
  $ lr(|S (z)|) < 1 , forall z < 0 upright(" and ") S (- oo) = lim_(z arrow.r - oo) S (z) = 0 $
]

Plugging $- oo$ int $S$ we obtain
$S (- oo) = 1 - bold(b)^top bold(frak(A))^(- 1) bold(1)$, which is equal to
zero if $bold(b)$ is equal to the last row of $bold(frak(A))$.

#theorem(number: "9.2.8.5", [Meta-theorem --- Convergence of fully discrete evolution])[
  Assume that 
  - the solution of the parabolic IBVP is "sufficiently smooth"
  - its spatial Galerkin finite element discretization relies on degree $p$ Lagrangian finite elements on uniformly shape-regular families of meshes
  - time stepping is based on an L-stable single step method of order $q$, 

  then we can expect an asymptotic behaviour
  of the total discretization error according to
  $ (tau sum_(j = 1)^M lr(|u (tau j) - u_h^((j))|)_(H^1 (Omega))^2)^(1 / 2) lt.eq C (h_(cal(M))^p + tau^q) $

]
Hence the total error is the spatial error plus the temporal error.

== Linear wave equations
<sub:linear-wave-equations>

In local form, the (linear) wave equation is given by
#neq($ rho (bx) frac(partial^2, partial t^2) u - div (sigma (bx) grad u) = f quad upright("in ") tilde(Omega) $) <eq:wave-strong>

Note the similarity to the heat equation @eq:heat-local. Since the wave equation is a second order ODE $accent(bold(u), dot.double) = bold(f (u))$,
two initial conditions are needed. In addition to the initial conditions @eq:heat-bc-ic, the initial velocity
$ frac(partial, partial t) u (bx , 0) = v_0 (bx) quad upright("for all ") bx in Omega $
is also needed.

We want to use the time stepping schemes we already know. To apply them,
the wave function can be converted into a first order ODE:
$ dot(u) & = v\
rho dot(v) & = div (sigma (bx) grad u) $ Remember
from Analysis that the particular wave equation
$frac(partial^2, partial t^2) u - c^2 frac(partial^2, partial x^2) u = 0$
in 1D results in the d’Alembert solution:
$ u (x , t) = 1 / 2 (u_0 (x + c t) + u_0 (x - c t)) + frac(1, 2 c) integral_(x - c t)^(x + c t) v_0 (s) dif s $
with $u_0$ and $v_0$ the initial conditions. Hence there is again the
concept of domain of dependence and domain of influence. This will be
important later. Furthermore, in the absence of a source term, as in the
simple case above, the solution will stay undamped. This corresponds to
#emph[conservation of total energy];.

We can formulate the variational problem:
#neq($ m (accent(u, dot.double) , v) + a (u , v) = 0 quad forall v in V_0 $) <eq:wave-variational>

#theorem(number: "9.3.2.16", "Energy conservation in wave propagation")[
  If $u$ solves @eq:wave-variational, then its energy is conserved, in the sense that
  $ 1 / 2 m (frac(partial, partial t) u , frac(partial, partial t) u) + 1 / 2 a (u , u) equiv upright("const") $
  where
  $1 / 2 m (frac(partial, partial t) u , frac(partial, partial t) u)$ can be understood as the 'kinetic' energy and $1 / 2 a (u , u)$ as the 'potential' energy.
]
#counter(heading).update((9,3,2))
=== Method of Lines

The method of lines gives rise to
$ bold(M) {frac(d^2, d t^2) arrow(mu) (t)} + bA arrow(mu) (t) & = arrow(phi) (t)\
arrow(mu) (0) = arrow(mu)_0 , frac(dif, dif t) arrow(mu) (0) & = arrow(nu)_0 $
Using $arrow(nu) = dot(arrow(mu))$, we can rewrite it to be a first
order ODE.
$ frac(dif, dif t) arrow(mu) & = arrow(nu)\
bold(M) frac(dif, dif t) arrow(nu) & = arrow(phi) (t) - bA arrow(mu)\
arrow(mu) (0) & = arrow(mu)_0 , arrow(nu) (0) = arrow(nu)_0 $ Remember
that in the case of $arrow(phi) equiv 0$, energy is conserved:
$ E_h (t) = 1 / 2 frac(dif, dif t) arrow(mu)^top bold(M) frac(dif, dif t) arrow(mu) + 1 / 2 arrow(mu)^top bA arrow(mu) equiv upright("const") $
So we would like a time stepping that preserves this. Such time stepping
schemes are called #emph[structure preserving];. One such timestepping
scheme is the #strong[Crank–Nicolson] one:
$ bold(M) frac(arrow(mu)^((j + 1)) - 2 arrow(mu)^((j)) + arrow(mu)^((j - 1)), tau^2) = - 1 / 2 bA (arrow(mu)^((j - 1)) + arrow(mu)^((j + 1))) + 1 / 2 (arrow(phi) (t_j - 1 / 2 tau) + arrow(phi) (t_j + 1 / 2 tau)) $
Another one would be the #strong[Störmer scheme];:
$ bold(M) frac(arrow(mu)^((j + 1)) - 2 arrow(mu)^((j)) + arrow(mu)^((j - 1)), tau^2) = - bA arrow(mu)^((j)) + arrow(phi) (t_j) $
For both of these #emph[second order] time stepping schemes, we need
$arrow(mu)^((- 1))$ to get $arrow(mu)^((1))$. Now the question is, where
do we get this from? It can be obtained with a special initial step,
using a symmetric (first order) difference quotient:
$ frac(dif, dif t) arrow(mu) (0) = arrow(nu)_0 arrow.r frac(mu^(\(1\)) - mu^((- 1)), 2 tau) = arrow(nu)_0 $
And finally there is the #strong[Leapfrog] timestepping. Using the
auxiliary variable
$arrow(nu)^((j + 1 \/ 2)) = frac(arrow(mu)^((j + 1)) - arrow(mu)^((j)), tau)$
and inserting this into the Störmer scheme results in
$ bold(M) frac(arrow(nu)^((j + 1)) - arrow(nu)^((j)), tau) & = - bA arrow(mu)^((j)) + arrow(phi) (t_j)\
frac(arrow(mu)^((j + 1)) - arrow(mu)^((j)), tau) & = arrow(nu)^((j + 1 \/ 2)) $
with the initial step
$arrow(nu)^((- 1 \/ 2)) + arrow(nu)^((1 \/ 2)) = 2 arrow(nu)_0$.
