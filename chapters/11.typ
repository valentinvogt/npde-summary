#import "../src/setup.typ": *
#show: thmrules

= Numerical Methods for Conservation Laws
<ch:conservation-laws>

#counter(heading).update((11, 1))
== Scalar Conservation Laws in 1D
<sub:scalar-conservation-laws-in-1d>

The goal of this chapter is to solve general conservation laws which are of the
form
#neq(
  $ frac(partial u, partial t) (x , t) + frac(partial, partial x) (f (u (x , t) , x)) = s (u (x , t) , x , t) . $,
) <eq:general-cons-law>
The flux $f : bb(R) times Omega arrow.r bb(R)$ can be a general function, which
can depend non-linearly on the solution $u$. Everything in this chapter will be
one-dimensional in space and time, so we have
$Omega subset.eq bb(R)$.

We usually consider the special case $s = 0$ and $f=f(u)$, known as the *Cauchy
problem*:
#neq(
  $ frac(partial u, partial t) + frac(partial, partial x) f (u) &= 0 \
  u(x,0)                                                      &= u_0 (x) $,
) <eq:cauchy-problem>
#mybox(
  "Example: Particle Model", ..unimportant,
)[
  We model cars as particles with position $x_i (t)$. Traffic speed is modeled by
  $ dot(x)_i (t) = v_"max" dot.op (1 - frac(Delta_0, Delta x_i (t))) , quad Delta x_i (t) = x_(i + 1) (t) - x_i (t) , $
  where $v_"max"$ is the maximum velocity and $Delta_0$ is the minimal distance
  between cars (i.e., the car length). Using some basic assumptions, we get a PDE:
  $ frac(partial u, partial t) (x , t) = frac(partial, partial x) u (1 - u) . $
  We define the car density $u(x)$ as the number of cars in an infinitesimal
  interval around $x$.
]

#counter(heading).step(level: 3)
=== Characteristics
<sub:characteristics>

We consider the Cauchy problem @eq:cauchy-problem. Then a characteristic curve
is defined as
#definition(
  number: "11.2.2.3", "Characteristic curve for 1D scalar conservation law",
)[
  $Gamma : [0 , T] arrow.r bb(R) times [0 , T]$ with $Gamma (tau) := (gamma (tau) , tau)$,
  such that $gamma$ satisfies
  $ frac(dif, dif t) gamma (tau) = f prime (u (gamma (tau) , tau)) $
  for $0 <= tau <= T$.
]

Generally, characteristic curves are lines along which information propagates.
This means $u (x , t)$ will only depend on the initial condition at $x_0$,
namely $u_0 (x_0)$, if there is a characteristic curve that starts in the point $x_0$ and
travels to the spacetime point $(x , t)$. One property of characteristic curves
is the following:

#lemma(
  number: "11.2.2.6", "Classical solution and characteristic curves",
)[
  Smooth solutions of @eq:cauchy-problem are constant along characteristic curves.
]

For example in the case of linear advection ($partial_t u + v partial_x u = 0$)
we can use this to solve the equation because $gamma (tau) = (x_0 + tau v)$,
which implies
$u (x , t) = u_0 (x - t v)$.

But this does not work if the solution is not smooth. For example, in the
traffic flow model above, the solution has a jump after a certain time and hence
this approach breaks down.

#counter(heading).step(level: 3)
=== Jump conditions and Riemann Problem
<sub:jump-conditions-and-riemann-problem>

The method of characteristics usually only works up to a certain point in time.
To get the solution for times after that, we first note that the solution will
usually have a discontinuity after the time where the method of characteristics
breaks down. So we study how the solution behaves at these jumps
(discontinuities).

The setting is as follows: We still study the Cauchy problem @eq:cauchy-problem.
We can derive that along jumps, the normal components must be continuous, which
leads to the
#definition(
  number: "11.2.4.2", [Rankine--Hugoniot (jump) condition],
)[
  $ dot(s) (u_l - u_r) = f (u_l) - f (u_r) $ where
  $dot(s) = frac(d gamma, d t)$ is the time derivative of the discontinuity curve $Gamma (t) = (gamma (t) , t) in bb(R) times [0 , T]$.
]

Note that this is useful because it allows us to compute the jump if we know $u_l$ and $u_r$.

The Riemann problem is given as
#definition(
  number: "11.2.5.1", "Riemann problem",
)[
  $ frac(partial u, partial t) + frac(partial f (u), partial x) = 0 $ and
  $ u_0 (x) = cases(
    delim: "{", u_l in bb(R) & upright("if ") x < 0, u_r in bb(R) & upright("if ") x > 0,

  ) $
  Note that $f$ can still be chosen to be any sufficiently smooth flux function.
]

Using the Rankine--Hugoniot jump condition we then get the following solution
for Riemann problems with a shock:
#lemma(
  number: "11.2.5.4", "Shock solution for Riemann problem",
)[
  For any two states $u_r , u_l in bb(R)$ the piecewise constant function
  $ u (x , t) := cases(
    delim: "{", u_l & upright("for ") x < dot(s) t, u_r & upright("for ") x > dot(s) t,

  ) & & quad quad dot(s) := frac(f (u_l) - f (u_r), u_l - u_r) , & & x in bb(R) , 0 < t < T $
  is a weak solution to the Riemann problem.
]

Note that the solution only holds if the equation implies a shock (jump). This
is the case if $f$ is convex and $u_l > u_r$ or if $f$ is concave and $u_r > u_l$.

If the jump only exists in the beginning we have a different solution
#lemma(
  number: "11.2.5.4", "Rarefaction solution for Riemann problem",
)[
  If $f in C^2 (bb(R))$ is strictly $cases(
    delim: "{", upright("convex and ") u_l < u_r, upright("concave and ") u_r < u_l,

  )$, then
  $ u (x , t) := cases(
    delim: "{", u_l & upright("for ") x < min { f prime (u_l) , f prime (u_r) } dot.op t, g (x / t) & upright("for ") min { f prime (u_l) , f prime (u_r) } < x / t < max { f prime (u_l) , f prime (u_r) }, u_r & upright("for ") x > max { f prime (u_l) , f prime (u_r) } dot.op t,

  ) $
  is a weak solution to the Riemann problem with $g := (f prime)^(- 1)$.
]

The question when to choose which of the two solution is answered by

#definition(
  number: "11.2.6.1", "Lax entropy condition",
)[
  Let $u$ be a weak solution of @eq:cauchy-problem, and a piecewise classical
  solution in neighborhood of $C^2$-curve $Gamma := (gamma (tau) , tau) , 0 <= tau <= T$,
  discontinuous across $Gamma$.

  $u$ satisfies the #emph[Lax entropy condition] in $(x_0 , t_0) in Gamma$ iff.
  $ f prime (u_l) > underbrace(frac(f (u_l) - f (u_r), u_l - u_r), dot(s)) > f prime (u_r) $
]

Now if $u$ satisfies the Lax entropy condition, then we have to pick the shock
solution. Otherwise we pick the rarefaction solution.

#counter(heading).update((11, 2, 6))
=== Properties of Entropy Solutions
<sub:properties-of-entropy-solutions>

The essential properties here are that with the propagation speed
$f prime (u)$, we find the #strong[domain of dependence] and the
#strong[domain of influence];, which is best illustrated by a picture and hence
we encourage the reader to look at the lecture document and the illustrations
below Theorem 11.2.7.3.

Moreover the second result is that the number of extrema of the solution is
non-increasing in time.

== Conservative Finite Volume (FV) Discretization
<sub:conservative-finite-volume>
=== Finite Difference Methods
<sub:finite-difference-methods>

Finite difference methods are probably the simplest methods for solving PDEs. We
just replace the spacial derivatives by some finite difference quotient for
example one of the following

#subtle-box[
  $ "Symmetric difference quotient" quad frac(partial f, partial x) &approx frac(f (x_0 + h) - f (x_0 - h), 2 h) \

  "Backward difference quotient" quad frac(partial f, partial x)  &approx frac(f (x_0) - f (x_0 - h), h) \

  "Forward difference quotient" quad frac(partial f, partial x)   &approx frac(f (x_0 + h) - f (x_0), h) $
]

Then we construct a solution by time stepping: given
$u (x , t_k)$ we compute $u (x , t_(k + 1))$ by using some Runge--Kutta
integrator.

With the example of finite difference methods, we observe that by the nature of
the problems we are studying in this chapter, the solutions have to be
constructed under consideration of the flux direction. That is, in the spirit of
characteristic curves we know that information propagates along curves in space
time. And if this curve advances, for example, from left to right in space, then
we need to use the forward difference quotient, because the backward difference
quotient will not contain the information of the flow direction.// TODO

=== Spatially Semi-Discrete Conservation Form
<sub:spatially-semi-discrete-conservation-form>

The method we will use in this section of the course is the Finite Volume
Method: we build a mesh by taking intervals around the spatial points in which
we approximate the solution $u$. And then we use the problem definition
$ frac(partial u, partial t) = - frac(partial, partial x) f (u) $ to derive in
the simplest case
#neq(
  $ frac(partial u (x_i , t_k), partial t)& approx frac(partial, partial t) 1 / h integral_(x_(j - 1 \/ 2))^(x_(j + 1 \/ 2)) u (x , t_k) dif x \
                                        & approx - 1 / h lr(
    ( F (u (x_j , t_k) , u (x_(j + 1) , t_k)) - F (u (x_(j - 1) , t_k) , u (x_j , t_k))), size: #150%,

  ) $,
) <eq:two-points-flux>

where $F (dot,dot)$ is an approximation of the flux function $f (u)$. Note that
in the above example $F$ only depends on two neighboring nodes, but this can
also be extended to use more nodes. Now we discretize this similarly to @ch:9,
where we let the coefficients $mu_i$ of the basis expansion of $u$ vary in time.
Here, we define
$ mu_j (t) &= u(x_j, t) \
bmu      &= (mu_1, dots, mu_N) $
With this, @eq:two-points-flux becomes (we omit the time argument of the $mu_i$)
#neq(
  $ frac(dif mu_j, dif t) (t) = - 1 / h (F (mu_j, mu_(j+1)) - F (mu_(j-1), mu_j)) $,
)<eq:two-points-flux-discrete>
We can simply plug the RHS of @eq:two-points-flux-discrete into a Runge--Kutta
method, so the only thing left to do is to find a suitable flux function $F$.

#definition(number: "11.3.3.5", "Consistent numerical flux function")[
  A numerical function $F : bb(R)^(m_l + m_r) arrow.r bb(R)$ is
  #emph[consistent] with the flux $f : bb(R) arrow.r bb(R)$ if
  $ F (u , u , dots.h , u) = f (u) , #h(2em) forall u in bb(R) $
]

Then we have the result that consistent numerical flux functions will
reconstruct the correct \"discrete shock speed\" when applied as explained
above. As this is such an important property, the finite volume method comes in
very handy to solve these Cauchy problems.

#counter(heading).step(level: 3)
=== Numerical Flux Functions
<sub:numerical-flux-functions>

This section now treats how to find suitable flux functions. There are several
options introduced starting with the simplest which basically corresponds to a
average of the two inputs in $F (u , w)$. But then it turns out that this flux
suffers from similar problems as the finite difference methods did.

One remedy for this is the #emph[Lax--Friedrichs / Rusanov] Flux which is useful
but flattens the edges of jumps (which is due to its construction with
additional diffusion).
#equation(
  number: "11.3.4.16", [Lax--Friedrichs / Rusanov Flux],
)[
  $ F_(L F) (v , w) = 1 / 2 (f (v) + f (w)) - 1 / 2 (w - v) max_(min { v, w } <= u <= max { v, w }) f prime (u) $
  The third term is the artificial diffusion term.
]
As pointed out before, the direction in which the information flows is crucial,
so an important idea to choose the right flux is to respect that. Moreover the
flux has to reproduce physical solution in the sense explained above when
studying two possible solutions for the Riemann problem.

The final flux we look at, which solves these problems, is the Godunov Flux.

#definition(
  number: "11.3.4.33", "Godunov Flux",
)[
  $ F_(G D) (v , w) = cases(
    delim: "{", min_(v lt.eq u lt.eq w) f (u) quad upright("if ") v < w, max_(w lt.eq u lt.eq v) f (u) quad upright("if ") v gt.eq w,

  ) $
]
#pagebreak(weak: true)
=== Monotone Schemes
<sub:monotone-schemes>

First, we define what it means for a flux to be monotone:
#definition(
  number: "11.3.5.5", "Monotonicity of flux functions",
)[
  A 2-point numerical flux function $F = F(v, w)$ is #emph[monotone] if

  $ - thick F "is increasing in its first argument, i.e."  &F(v, w) <= F(v+Delta v, w) quad &forall w &in bb(R) \
  - thick F "is decreasing in its second argument, i.e." &F(v, w) <= F(v, w+Delta w) quad &forall v &in bb(R)$
]<def:monotone-flux>
In @sub:properties-of-entropy-solutions, we mentioned that the number of extrema
of an analytical solution does not increase over time. Here, we show that the
two main fluxes (LF and Godunov) we derived in the previous chapter are monotone
and therefore both have this property. This is established by the following two
lemmas:
#lemma(
  number: "11.3.5.8", [Monotonicity of Lax--Friedrichs and Godunov flux],
)[
  For any continuously differentiable flux function $f$ the associated
  Lax--Friedrichs flux and Godunov flux are monotone.
]
#v(-1cm)
#lemma(
  number: "11.3.5.13", "Non-oscillatory monotone semi-discrete evolutions",
)[
  If $bmu = bmu (t)$ solves the two-point flux equation
  @eq:two-points-flux-discrete with a monotone numerical flux and $bmu (0)$ has
  finitely many local extrema, then the number of local extrema of $bmu (t)$ cannot
  be larger than that of $bmu (0)$.
]

== Timestepping for Finite-Volume Methods
<sub:timestepping-for-fv>

As explained before, once we have chosen the numerical Flux, we just need to
apply Runge--Kutta integration. This subsection studies some conditions that
have to be considered when applying Runge--Kutta -- in particular, constraints
on choosing the timestep size $tau$.

#definition(
  number: "11.4.2.5", "Numerical domain of dependence",
)[
  Consider the explicit fully discrete evolution
  $bmu^((k + 1)) := fvH (bmu^((k)))$ on a uniform spatio-temporal mesh ($x_j = h j , j in bb(Z) , t_k = k tau , k in bb(N)_0$)
  with
  #neq(
    $ exists m in bb(N)_0 : (fvH (bmu))_j = fvH (mu_(j - m) , dots.h , mu_(j + m)) , j in bb(Z) . $,
  ) <eq:fv-stencil-width>
  Then the *numerical domain of dependence* is given by
  $ D_h^(-) (x_j , t_k) := { (x_n , t_l) in bb(R) times [0 , t_k] : j - m (k - l) lt.eq n lt.eq j + m (k - l) } $
]
$fvH$ is an operator which gives $bmu$ at the next timestep. It incorporates
both the numerical flux and the Runge--Kutta method. If we have a first-order
time integrator and the numerical flux function is $F=F(mu_(j-m), dots.h, mu_(j+m))$,
then $fvH = fvH(mu_(j-m), dots.h, mu_(j+m))$.

#subtle-box[
  *Example:* Let's say we are solving @eq:two-points-flux-discrete and call the
  RHS $R(mu_j)$
  $ frac(dif mu_j, dif t) (t) = - 1 / h (F (mu_j, mu_(j+1)) - F (mu_(j-1), mu_j)) = R(mu_j) $
  and we apply explicit Euler:
  $ (fvH (bmu))_j &= mu_j + tau R(mu_j)\
                &= mu_j - tau / h (F (mu_j, mu_(j+1)) - F (mu_(j-1), mu_j)) $
  If we compare this to @eq:fv-stencil-width, we see that $m = 1$ in this case.
]
// TODO Elaborate

The following kind of condition appears over and over in numerical integration
and gives an upper bound for the timestep size $tau$: If this upper bound is
respected, the numerical solution is stable.
#definition(
  number: "11.4.2.11", [Courant--Friedrichs--Lewy (CFL) condition],
)[
  An explicit local fully discrete evolution
  $bmu^((t + 1)) := fvH (bmu)$ on uniform spatio-temporal mesh ($x_j = h j , j in bb(Z) , t_k = k tau , k in bb(N)$)
  satisfies the CFL condition, if the convex hull of its numerical domain of
  dependence contains the maximal analytical domain of dependence
  $ D^(-) (x_j , t_k) subset "convex"thin (D_h^(-) (x_j , t_k)) quad forall j , k . $
]
Applied to our problem, this means
#equation(
  number: "11.4.2.12", "CFL condition for explicit fully discrete evolution",
)[
  $ tau / h lt.eq frac(m, max { lr(|dot(s)_min|) , lr(|dot(s)_max|) }) . $
]

#counter(heading).update((11, 4, 3))
=== Convergence of Fully Discrete FV Method
<sub:convergence-of-fully-discrete-fv>

The *consistency error* is defined as follows:
$ eps := max_j {F(u(x_j, t), u(x_(j+1), t)) - f(u(x_(j+1/2), t))}, $
if we assume $u$ to be an exact solution of If $eps = Order(h^q)$, the flux is
called _q-th order consistent_. $q$ is then the order of convergence of the FV
method.

We get at most order one convergence in the maximum and the
$L^1$ norm. This can be seen by the following fact:
#mybox("Order barrier for monotone numerical fluxes")[
  Monotone numerical fluxes (@def:monotone-flux) are at most first order
  consistent.
]
#pagebreak(weak: true)
== Higher-Order FV Schemes
<sub:higher-order-fv-schemes>
For the FV Method to achieve a convergence of order $ > 1$, we need fluxes which
are consistent of order $ > 1$. Instead of coming up with completely new fluxes,
we use the ones from @sub:numerical-flux-functions and only modify the arguments
given to them.

We choose the approach to go from piecewise-constant to *discontinuous
piecewise-linear* solutions. Instead of approximating $f(u(x_(j+1/2)))$ by $F(mu_(j), mu_(j+1))$,
we use $F(nu^+_(j),nu^-_(j+1))$.

#grid(
  columns: (0.75fr, 0.2fr), column-gutter: 1em, [
    $nu^+_j,nu^-_j$ are the values of the linear approximation of $u$ at the
    boundaries of cell $j$ (see image). The discrete evolution equation becomes
    #neq(
      $ frac(d mu_j, d t) (t) = - 1 / h (F (nu_j^+ (t), nu_(j+1)^- (t)) - F (nu_(j-1)^+ (t), nu_j^- (t)) ) $,
    )<eq:linear-reconstruction-evolution>
  ], [
    #image("../images/linear_reconstruction.png", width: 80%)
  ],
)
where $nu_j^+,nu_j^-$ are the values of the linear approximation of $u$ at the
cell boundaries. The linear approximation on cell $j$ is defined by the slope $sigma_j$.
So this method boils down to finding formulas for the slopes:
$ sigma_j = "slope"(mu_(j-1), mu_j, mu_(j+1)) $

#counter(heading).step(level: 3)
=== Slope limiting
Finding a linear approximation is formally expressed as a *reconstruction
operator*:
$ recop = "discontinuous pw. linear approximation of " u $
Note that this operator as a mapping between vector spaces is not linear.
However, we call it a *linear reconstruction* because the result is a
(piecewise) linear function. We get the following formulas which we can plug
into @eq:linear-reconstruction-evolution:
$ v^+_(j) = recop (x^-_(j+1/2)) " and " v^-_j = recop (x^+_(j-1/2)) $
where $x^+$ denotes the limit from the right and $x^-$ from the left. Note that
we always assume $recop (x_j) = mu_j$, i.e., the linear reconstruction has the
same values as the piecewise constant function at the cell centers.

Given slopes $sigma_j$, we can now define the linear reconstruction operator:
#v(0.2cm)
#subtle-box(
  width: 70%,
)[
  #v(-0.2cm)
  $ recop (x) = mu_j + sigma_j (x - x_j) upright(" for ") x in openint(x_(j-1/2), x_(j+1/2)) $
  #v(-0.2cm)
]

We introduce an important property of the reconstruction operator:
#definition(
  number: "11.5.2.1", "Monotonicity-preserving linear reconstruction",
)[
  A linear reconstruction operator $upright(R)_M$ is called #emph[monotonicity-preserving] if
  $ mu_j <= mu_(j+1) &==> recop upright("is non-decreasing") in openint(x_j, x_(j+1)) \
  mu_j >= mu_(j+1) &==> recop upright("is non-increasing") in openint(x_j, x_(j+1)) $
]

This ensures that no new extrema are introduced by the numerical scheme.

Finally, we present a reconstruction operator which is monotonicity-preserving:
#definition(
  number: "11.5.2.9", "Minmod reconstruction",
)[

  $ sigma_j (mu_(j-1), mu_j, mu_(j+1)) = "minmod" (frac(mu_(j+1) - mu_j, x_(j+1) - x_j), frac(mu_j - mu_(j-1), x_j - x_(j-1))) $
  where $"minmod"$ is defined as
  $ "minmod" (a, b) = cases(
    delim: "{", a \, quad & upright("if ") quad a b > 0\, abs(a) &< abs(b) \
    b \,      & upright("if ") quad a b < 0\, abs(b) &< abs(a) \
    0 \,      & upright("otherwise"),

  ) $
]
