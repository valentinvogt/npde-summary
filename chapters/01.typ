#import "../src/setup.typ": *
#show: thmrules

= Second-Order Scalar Elliptic Boundary Value Problems
<ch:second-order-scalar-elliptic-boundary-value-problems>
#counter(heading).step(level: 2)
== Quadratic Minimization Problems
<sub:quadratic-minimization-problems>
#v(0.2cm)
In the following, let $V$ be a Vector Space over $bb(R)$. \
#mybox("Linear forms")[
  $ell : V arrow.r bb(R)$ is a #emph[linear form / linear functional] $arrow.l.r.double.long$
#neq($ ell (alpha u + beta v) = alpha ell (u) + beta ell (v) wide forall u , v in V , forall alpha , beta in RR $)
]
#v(-1cm)

#let bs = colMath("s", rgb("2900d3"))
#let bt = colMath("t", rgb("2900d3"))
#mybox("Bilinear forms")[
  $a : V times V arrow.r bb(R)$ is a #emph[bilinear form] $arrow.l.r.double.long$
#neq($ a \( &bs u_1 + u_2 , bt v_1 + v_2 \)\
 & = bs bt dot.op a (u_1 , v_1) + bs a (u_1 , v_2) + bt a (u_2 , v_1) + a (u_2 , v_2) forall u_1 , u_2 , v_1 , v_2 in V , forall bs , bt in bb(R) $)
]
#v(-1cm)
#mybox("Positive definiteness")[
  A bilinear form $a : V times V arrow.r bb(R)$ is #emph[positive definite] if
    $ u in V \\ { bold(0) } arrow.l.r.double.long a (u , u) > 0 $ 
  It is #emph[positive semi-definite] if
    $ a (u , u) & gt.eq 0 quad forall u in V $
]
#v(-1cm)

#mybox("Quadratic Functional")[
  A #emph[quadratic functional] $J : V arrow.r bb(R)$ is defined by
  #neq($ J (u) := 1 / 2 a (u , u) - l (u) + c, quad u in V $) 
  where $a : V times V arrow.r bb(R)$ a symmetric bilinear form, $ell : V arrow.r bb(R)$ a linear form and $c in bb(R)$. 
]
#v(-1cm)

#mybox("Continuity of linear form")[
  A linear form $ell : V arrow.r bb(R)$ is #emph[continuous / bounded] on $V$, if
  #neq($ exists C > 0 quad lr(|ell (v)|) lt.eq C norm(v) quad forall v in V $)
]

== Sobolev Spaces
<sub:sobolev-spaces>
When we solve a quadratic minimization problem, i.e., a quadratic
function for which we search a minimizer, we first need to define the
space of functions in which we want to look for the solution. For
example, in Physics, we generally want the solution (function that
describes e.g. the shape of an elastic string) to be continuous. So it
makes no sense to look for a minimizer with jumps. \
\
To formulate this mathematically we need the #strong[Sobolev space];.
For functions in a Sobolev space, it is ensured that the bilinear form
in the quadratic functional is well defined (i.e. finite). Hence the
mathematical space in which we look for minimizers is determined by the
given quadratic functional. To select the space for your problem, follow
the guideline ... #emph[Choose the largest space such that the problem
is well defined];.

#mybox("Sobolev Spaces")[
  $H^1_0 (Omega)$ is a vector space with norm 
  $ |v|_(H^1) := (integral_Omega norm(grad v)^2 dif bx)^(1 / 2) $
  $H^1 (Omega)$ is another vector space with norm 
  $ norm(v)^2_(H^1) := norm(v)^2_(L^2) + |v|^2_(H^1) $
  Note that $|dot|_(H^1)$ is not a norm on the space $H^1 (Omega)$, but a seminorm. Both spaces contain all functions for which the norm is finite (and, in the case of $H^1_0 (Omega)$, which are 0 on $partial Omega$).
]
*Alternative notation* for norms includes $norm(dot)_0$ for $norm(dot)_(L^2)$ and $|dot|_1$ for $|dot|_(H^1)$.

If the quadratic minimization problem is well defined, we get the
following lemma for existence and uniqueness of minimizers:

#theorem(number: "1.3.3.6", "Existence of minimizers in Hilbert spaces")[
  On a real Hilbert space $V$ with norm $norm(dot)_a$ for any
  $norm(dot)_a$-bounded linear functional
$ell : V arrow.r bb(R)$ the quadratic minimization problem
#neq($ u_(\*) & = op("argmin", limits: #true)_(v in V) J (v)\
J (v) & := 1 / 2 norm(v)_a^2 - ell (v) $)
 has a unique
  solution.
]

Note that here, we use the bilinear form to define the norm
$a (u , u) = norm(u)_a^2$. The main point is that in
(energy) minimization problems, the bilinear form of the quadratic
minimization problem can be seen as the norm of some Sobolev space. This
then leads to a solution if we check boundedness of the linear form.

For checking boundedness we can often use Cauchy-Schwarz and
Poincaré-Friedrichs inequalities.

== Linear Variational Problem
<sub:linear-variational-problem>
#definition(number: "1.4.1.6", "Linear Variational Problem")[
  Let $V$ be a vector (function) space, $mhat(V) subset V$ an affine space, and $V_0 subset V$ the associated subspace. The equation
#neq($ #text("Find") u in mhat(V) med #text("such that") a (u , v) = ell (v) quad forall v in V_0 $)
is called a (generalized) #emph[linear variational problem];, if

- $a : V times V_0 arrow.r bb(R)$ is a bilinear form

- $ell : V_0 arrow.r bb(R)$ is a linear form
]
Knowing that a solution exists is of course not enough. And solving a
minimization problem over infinite-dimensional spaces is not an easy
task. So we reformulate the problems in a linear variational form which
is then already pretty close to what we can solve numerically. For this,
we have the following equivalence:

#theorem(number: "1.4.1.8", "Equivalence of quadratic
minimization problem and linear variational problem")[
  For a (generalized) quadratic functional $J (v) = 1 / 2 a (v , v) - ell (v) + c$ on a vector space $V$ and with a symmetric positive definite bilinear form $a : V times V arrow.r bb(R)$ the following is equivalent: 

  - The quadratic minimization problem for $J (v)$ has the unique minimizer $u_(\*) in mhat(V)$ over the affine subspace $mhat(V) = g + V_0 , g in V$
  - The linear variational problem $ u in mhat(V) quad a (u , v) = ell (v) &quad forall v in V_0 $ has the unique solution $u_(\*) in mhat(V).$
]<thm:variational-problem-equiv>
Note that the test space $V_0$ and the trial space $mhat(V)$ can be
different. When our problem has Dirichlet boundary conditions, we
incorporate them into the trial space $mhat(V)$, which is also the space
from which we pick a solution. In the test space $V_0$ we do not need to
respect the boundary conditions. Instead, we set test functions to zero
where the boundary data is given (\"Don’t test where the solution is
known!\").

== Boundary Value Problems
<sub:boundary-value-problems>
#lemma(number: "1.5.2.1", "General product rule")[
  For all $bold(j) in (C^1 (overline(Omega)))^d , v in C^1 (overline(Omega))$ holds
  #neq($ div (bold(j) v) = v div bold(j) + bold(j) dot.op grad v quad upright("in") Omega $)
] <thm:general-product-rule>

#lemma(number: "1.5.2.4", "Gauss' Theorem")[
  Let $bold(n) : partial Omega arrow.r bb(R)^d$ denote the exterior unit normal vector field on $partial Omega$ and $d S$ denote integration over a surface. We have
  #neq($ integral_Omega div thin bold(j (x)) dif bx = integral_(partial Omega) bold(j (x) dot.op n (x)) dif S (x) quad forall bold(j) in (C_(upright(p w))^1 (overline(Omega)))^d $)
] <thm:gauss-theorem>

#lemma(number: "1.5.2.7", "Green's first formula")[
  For all vector fields $bold(j) in (C^1 (overline(Omega)))^d$ and functions $v in C^1 (overline(Omega))$ holds
  #neq($ integral_Omega bold(j) dot.op grad v dif bx = - integral_Omega div thin bold(j) thin v dif bx + integral_(partial Omega) bold(j dot.op n) v thin dif S $)
] <thm:greens-formula>

#lemma(number: "1.5.3.4", "Fundamental lemma of the calculus of variations")[
  Let $f in L^2 (Omega)$ satisfy
  #neq($ integral_Omega f (bx) v (bx) dif bx = 0 quad forall v in C_0^oo (Omega), $)
  then $f equiv 0$.
] <thm:variational-calculus-fundamental-lemma>

We have seen equivalence of minimization problem of a quadratic
functional and linear variational problem. The variational problem is called the #strong[weak form];, we can transform it (with extra smoothness
requirements) into the problem's #strong[strong form];, i.e., into an
elliptic BVP, mainly with the help of the above lemmas.

#counter(heading).step(level: 2)

== Boundary Conditions
<sub:boundary-conditions>
For 2nd-order elliptic BVPs we need boundary conditions to get a unique
solution. To be more precise, we need #strong[exactly one] of the
following boundary conditions on every part of $partial Omega$

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

== Second-Order Elliptic Variational Problems
<sub:second-order-elliptic-variational-problems>
We have seen how we can get from a minimization problem via a
variational problem to a BVP. Now we want to move in the opposite
direction, from a PDE and its boundary conditions we want to get to a
variational problem. This can again be done using the lemmas from
section 1.5 and considering the boundary conditions to choose a suitable
(Sobolev) function space.

For Neumann problems there is a #strong[compatibility condition];. If we
choose test function $v equiv 1$ we get the requirement
$ - integral_(partial Omega) h thin dif S = integral_Omega f thin dif bx $
for the existence of solutions. Additionally, the solution of Neumann
problems is unique only up to constants. To address this we can use the
constrained function space
$ H_(\*)^1 (Omega) := { v in H^1 (Omega) : integral_Omega v thin dif bx = 0 } $

#theorem(number: "1.8.0.20", title: "Theorem", "Second Poincaré-Friedrichs inequality")[
  If $Omega subset bb(R)^d$ is bounded and connected, then
  #neq($ exists C = C (Omega) > 0 : norm(u)_0 lt.eq C "diam" #h(-0.1pt) (Omega) thin norm(grad u)_0 quad forall u in H_(\*)^1 (Omega) $)
] <thm:poincare-friedrichs>
This theorem tells us that (under some conditions), the $L^2$ norm of functions from this space is bounded by the $H^1$-seminorm.

== Essential and Natural boundary Conditions
<sub:essential-and-natural-boundary-conditions>
Essential boundary conditions are boundary conditions which have been
imposed directly on the trial space, i.e. Dirichlet BC. On the other
hand, Neumann BC are only enforced through the variational equation. They are callednatural boundary conditions.

- #strong[Admissible Dirichlet Data];: Dirichlet boundary values need to
  be continuous.

- #strong[Admissible Neumann Data];: $h$ needs to be in $L^2 (Omega)$
  (can be discontinuous)

#theorem(number: "1.9.0.19", title: "Theorem", "Multiplicative trace inequality")[
  #neq($ exists C = C (Omega) > 0 : norm(u)_(L^2(partial Omega)) lt.eq C norm(u)_(L^2(Omega)) dot.op norm(u)_(H^1(Omega)) quad forall u in H^1 (Omega) $)
] <thm:mult-trace-inequality>
This theorem is frequently needed when dealing with integrals over the boundary.
