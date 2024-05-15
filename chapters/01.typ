#import "../src/setup.typ": *
#show: thmrules

= Second-Order Scalar Elliptic Boundary Value Problems
<ch:01>
#counter(heading).step(level: 2)
== Quadratic Minimization Problems
<sub:quadratic-minimization-problems>
#v(0.2cm)

In the following, let $V$ be a vector space over $bb(R)$.\
#mybox(
  "Linear forms",
)[
  $ell : V arrow.r bb(R)$ is a #emph[linear form / linear functional] $arrow.l.r.double.long$
  #neq(
    $ ell (alpha u + beta v) = alpha ell (u) + beta ell (v) wide forall u , v in V , forall alpha , beta in RR $,
  )
]
#v(-1cm)

#let bs = colMath("s", rgb("2900d3"))
#let bt = colMath("t", rgb("2900d3"))
#mybox(
  "Bilinear forms",
)[
  $a : V times V arrow.r bb(R)$ is a #emph[bilinear form] $arrow.l.r.double.long$
  #neq(
    $ a \( &bs u_1 + u_2 , bt v_1 + v_2 \)\
         & = bs bt dot.op a (u_1 , v_1) + bs a (u_1 , v_2) + bt a (u_2 , v_1) + a (u_2 , v_2) forall u_1 , u_2 , v_1 , v_2 in V , forall bs , bt in bb(R) $,
  )
]
#v(-1cm)
#mybox(
  "Positive definiteness",
)[
  A bilinear form $a : V times V arrow.r bb(R)$ is #emph[positive definite] if
  $ u in V \\ { bold(0) } arrow.l.r.double.long a (u , u) > 0 $
  It is #emph[positive semi-definite] if
  $ a (u , u) & gt.eq 0 quad forall u in V $
]
#v(-1.1cm)

#mybox(
  "Quadratic Functional",
)[
  A #emph[quadratic functional] $J : V arrow.r bb(R)$ is defined by
  #neq($ J (u) := 1 / 2 a (u , u) - ell (u) + c, quad u in V $)
  where $a : V times V arrow.r bb(R)$ a symmetric bilinear form, $ell : V arrow.r bb(R)$ a
  linear form and $c in bb(R)$.
]
#v(-1cm)

#mybox(
  "Continuity of linear form",
)[
  A linear form $ell : V arrow.r bb(R)$ is #emph[continuous / bounded] on $V$, if
  #neq(
    $ exists C > 0 quad lr(|ell (v)|) lt.eq C norm(v) quad forall v in V, $,
  ) <eq:continuity-linear-form>
  where $norm(dot)$ is a norm on $V$.
]

== Sobolev Spaces
<sub:sobolev-spaces>

When we solve a minimization problem, we first need to define the space of
functions in which we look for the solution. For example, in Physics, we
generally want the solution to be continuous. E.g, a function describing the
shape of an elastic string should not have jumps.

It turns out that the correct space to describe our minimization problems is the #strong[Sobolev space];.
For functions $u$ in a Sobolev space, the bilinear form $a$ in the quadratic
functional is well defined (i.e., $a(u,u)<oo$). Hence the space in which we look
for minimizers is determined by the given quadratic functional. To select the
space for your problem, follow the guideline ...
#emph[Choose the largest space such that the problem is well defined];.

#mybox(
  "Sobolev Spaces",
)[
  $H^1_0 (Omega)$ is a vector space with norm
  $ |v|_(H^1) := (integral_Omega norm(grad v)^2 dif bx)^(1 / 2) $
  $H^1 (Omega)$ is another vector space with norm
  $ norm(v)^2_(H^1) := norm(v)^2_(L^2) + |v|^2_(H^1) $
  Note that $|dot|_(H^1)$ is not a norm on the space $H^1 (Omega)$, but a
  seminorm.

  Both spaces contain all functions for which the norm is finite (and, in the case
  of $H^1_0 (Omega)$, which are 0 on $partial Omega$).
]
*Alternative notation* for norms includes $norm(dot)_0$ for $norm(dot)_(L^2)$ and $|dot|_1$ for $|dot|_(H^1)$.

If the quadratic minimization problem is well defined, we get the following
lemma for existence and uniqueness of minimizers:

#theorem(
  number: "1.3.3.6", "Existence of minimizers in Hilbert spaces",
)[
  On a real Hilbert space $V$ with norm $norm(dot)_a$ for any
  $norm(dot)_a$-bounded linear functional $ell : V arrow.r bb(R)$, the quadratic
  minimization problem
  #neq($ u_(\*) & = op("argmin", limits: #true)_(v in V) J (v)\
  J (v)  & := 1 / 2 norm(v)_a^2 - ell (v) $)
  has a unique solution.
] <thm:existence-minimizer-hilbert>

Note that here, we use the bilinear form to define the norm
$norm(u)_a = sqrt(a (u , u))$. The main point is that we can see the bilinear
form of the quadratic minimization problem as the norm of some Sobolev space.
The above theorem guarantees that a solution exists in this space if the linear
form is bounded.

For checking boundedness we can often use Cauchy--Schwarz
(@eq:cauchy-schwarz-integrals) and Poincaré--Friedrichs
(@thm:poincare-friedrichs).

#pagebreak(weak: true)
== Linear Variational Problem
<sub:linear-variational-problem>

#definition(
  number: "1.4.1.6", "Linear Variational Problem",
)[
  Let $V$ be a vector (function) space, $mhat(V) subset V$ an affine space, and $V_0 subset V$ the
  associated subspace. The equation
  #neq(
    $ #text("Find") u in mhat(V) med #text("such that") a (u , v) = ell (v) quad forall v in V_0 $,
  ) <eq:linear-variational-problem>
  is called a (generalized) #emph[linear variational problem];, if

  - $a : V times V_0 arrow.r bb(R)$ is a bilinear form

  - $ell : V_0 arrow.r bb(R)$ is a linear form
]

@thm:existence-minimizer-hilbert tells us that the minimization problem has a
solution, but knowing that a solution exists is of course not enough: We want to
find it, but an infinite-dimensional minimization problem is hard to solve. To
make it easier, we reformulate the problems in a linear variational form
@eq:linear-variational-problem, which is quite close to something we can solve
numerically. To do this transformation, we use the following equivalence:

#theorem(
  number: "1.4.1.8", "Equivalence of quadratic
            minimization problem and linear variational problem",
)[
  For a (generalized) quadratic functional $J (v) = 1 / 2 a (v , v) - ell (v) + c$ on
  a vector space $V$ and with a symmetric positive definite bilinear form $a : V times V arrow.r bb(R)$ the
  following is equivalent:

  - The quadratic minimization problem for $J (v)$ has the unique minimizer $u_(\*) in mhat(V)$ over
    the affine subspace $mhat(V) = g + V_0 , g in V$.
  - The linear variational problem $ u in mhat(V) quad a (u , v) = ell (v) &quad forall v in V_0 $ has
    the unique solution $u_(\*) in mhat(V).$
]<thm:variational-problem-equiv>

Note that the trial space $mhat(V)$, from which we pick a solution, and the test
space $V_0$ can be different. For an example of different trial and test spaces,
see @sub:boundary-conditions.

#pagebreak(weak: true)
== Boundary Value Problems
<sub:boundary-value-problems>

#lemma(
  number: "1.5.2.1", "General product rule", ..unimportant,
)[
  For all $bold(j) in (C^1 (overline(Omega)))^d , v in C^1 (overline(Omega))$ holds
  #neq(
    $ div (bold(j) v) = v div bold(j) + bold(j) dot.op grad v quad upright("in") Omega $,
  )
] <thm:general-product-rule>

#lemma(
  number: "1.5.2.4", "Gauss' Theorem",
)[
  Let $bold(n) : partial Omega arrow.r bb(R)^d$ denote the exterior unit normal
  vector field on $partial Omega$ and $d S$ denote integration over a surface. We
  have
  #neq(
    $ integral_Omega div bold(j (x)) dif bx = integral_(partial Omega) bold(j (x) dot.op n (x)) dif S (bx) quad forall bold(j) in (C_(upright(p w))^1 (overline(Omega)))^d $,
  )
] <thm:gauss-theorem>

#lemma(
  number: "1.5.2.7", "Green's first formula",
)[
  For all vector fields $bold(j) in (C^1_"pw" (overline(Omega)))^d$ and functions $v in C^1_"pw" (overline(Omega))$ holds
  #neq(
    $ integral_Omega bold(j) dot.op grad v dif bx = - integral_Omega div bold(j) thin v dif bx + integral_(partial Omega) bold(j dot.op n) thin v dif S $,
  )
] <thm:greens-formula>

#lemma(
  number: "1.5.3.4", "Fundamental lemma of the calculus of variations",
)[
  Let $f in L^2 (Omega)$ satisfy
  #neq(
    $ integral_Omega f (bx) v (bx) dif bx = 0 quad forall v in C_0^oo (Omega), $,
  )
  then $f equiv 0$.
] <thm:fund-lemma>

We have seen that minimizing a quadratic functional is equivalent to solving a
linear variational problem @eq:linear-variational-problem. The variational
problem is called the #strong[weak form];. We can transform it (with extra
smoothness requirements) into the problem's #strong[strong form];, an elliptic
BVP (PDE with boundary conditions).
#tip-box(
  "Weak to strong",
)[
  + Use @thm:greens-formula to get rid of derivatives on $v$ (e.g. turn $grad u dot grad v$ into $-div(grad u) +...$
  + Use properties of the test space (usually that $v=0$ on $partial Omega$) to get
    rid of boundary terms
  + Use @thm:fund-lemma to remove the integrals and test functions
]

#pagebreak(weak: true)
#counter(heading).step(level: 2)
== Boundary Conditions
<sub:boundary-conditions>

For 2nd-order elliptic BVPs we need boundary conditions to get a unique
solution. To be more precise, we need #strong[exactly one] of the following
boundary conditions on every part of $partial Omega$

#mybox("Main boundary conditions for 2nd-order elliptic BVPs")[
  - #strong[Dirichlet]: $u$ is fixed to be
    $g : partial Omega arrow.r bb(R)$
    $ u = g quad upright("on") thin partial Omega $

  - #strong[Neumann]: the flux $bold(j) = - kappa (bx) grad u$
    through $partial Omega$ is fixed with
    $h : partial Omega arrow.r bb(R)$
    $ bold(j dot.op n) = - h quad upright("on") thin partial Omega $

  - #strong[Radiation]: flux depends on $u$ with an increasing function
    $Psi : bb(R) arrow.r bb(R)$
    $ bold(j dot.op n) = Psi (u) quad upright("on") thin partial Omega $
]

In the weak form, Dirichlet conditions have to be imposed directly on the
*trial* space. The test space needs to be set to 0 wherever Dirichlet conditions
are given ("_Don't test where the solution is known_"). For example, trial and
test spaces for a standard Dirichlet problem are
$ V   &= Set(u in C^1(Omega)&,   && u=g "on" partial Omega) \
V_0 &= Set(v in C^1(Omega)&,   && v=0 "on" partial Omega) $

Dirichlet BCs are called #strong[essential boundary conditions];.

Neumann conditions, which are only enforced through some term in the variational
equation, are called #strong[natural boundary conditions];.

There are some constraints on the boundary data:
#subtle-box[
  - #strong[Admissible Dirichlet Data];: Dirichlet boundary values need to be
    continuous.

  - #strong[Admissible Neumann Data];: $h$ needs to be in $L^2 (Omega)$
    (can be discontinuous)
]

The following theorem is frequently needed when dealing with integrals over the
boundary:
#theorem(
  number: "1.9.0.19", title: "Theorem", "Multiplicative trace inequality",
)[
  #neq(
    $ exists C = C (Omega) > 0 : norm(u)_(L^2(partial Omega)) lt.eq C norm(u)_(L^2(Omega)) dot.op norm(u)_(H^1(Omega)) quad forall u in H^1 (Omega) $,
  )
] <thm:mult-trace-inequality>

#pagebreak(weak: true)
== Second-Order Elliptic Variational Problems
<sub:second-order-elliptic-variational-problems>

We have seen how we can get from a minimization problem via a variational
problem to a BVP. Now we want to move in the opposite direction: from a PDE with
boundary conditions to a variational problem.
#tip-box(
  "Strong to weak",
)[
  + Test the PDE with (multiply by $v$) and integrate over $Omega$
  + Use @thm:greens-formula to "shift" one derivative from $u$ to $v$ (e.g., from $-div(grad u)$ to $grad u dot grad v + ...$)
  + Use Neumann BC on boundary terms ($grad u dot n = h$)
  + Pick Sobolev trial/test spaces $V,V_0$ such that

    - $a(u,u)$ is finite for $u in V,V_0$
    - boundary conditions are satisfied ($u=g$ in $V$ $=>$ $v=0$ in $V_0$)

    To fulfill the first condition, we can define the "base" space for both trial
    and test as $Set(v, a(v,v)<oo)$, which is equal to $H^1$ for the usual $Delta u = f$ problem.
    If there are extra (e.g., boundary) terms in $a$, try to bound these with the $H^1$ norm.
]
For Neumann problems there is a #strong[compatibility condition];. If we choose
test function $v equiv 1$ we get the requirement
$ - integral_(partial Omega) h dif S = integral_Omega f dif bx $
for the existence of solutions. Additionally, the solution of Neumann problems
is unique only up to constants. To address this we can use the constrained
function space
$ H_(\*)^1 (Omega) := { v in H^1 (Omega) : integral_Omega v dif bx = 0 } $

#theorem(
  number: "1.8.0.20", title: "Theorem", [Second Poincaré--Friedrichs inequality],
)[
  If $Omega subset bb(R)^d$ is bounded and connected, then
  #neq(
    $ exists C = C (Omega) > 0 : norm(u)_0 lt.eq C "diam" #h(-0.1pt) (Omega) thin norm(grad u)_0 quad forall u in H_(\*)^1 (Omega) $,
  )
] <thm:poincare-friedrichs>
This theorem tells us that (under some conditions), the $L^2$ norm of functions
from this space is bounded by the $H^1$-seminorm.