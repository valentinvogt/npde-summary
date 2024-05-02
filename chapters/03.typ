#import "../src/setup.typ": *
#import "../src/theorems.typ": *
#show: thmrules

= FEM: Convergence and Accuracy
<ch:convergence-and-accuracy>
== Abstract Galerkin Error Estimates
<sub:abstract-galerkin-error-estimates>

The main takeaway here is that the solution given by the Galerkin method is the
best one (with respect to the energy norm) in the chosen discrete subspace. This
is formalized by the following lemma:

#theorem(
  number: "3.1.3.7", "Cea's Lemma",
)[
  Under some assumptions that guarantee the existence of a unique solution we have
  #neq($ norm(u - u_h)_a = inf_(v_h in V_(0 , h)) norm(u - v_h)_a, $)
  where $u$ is the exact solution and $u_h$ is the Galerkin solution.
]<thm:ceas-lemma>

Next we want to discuss the types of refinement, i.e., the steps we can take to
increase the accuracy of our method.

#mybox(
  "h-Refinement",
)[
  Replace the mesh $cal(M)$ (underlying $V_(0 , h)$) with a finer mesh $cal(M) prime$ (underlying
  larger discrete trial space $V prime_(0 , N prime)$).
] <concept:refinement>
#v(-1cm)
#mybox(
  "p-Refinement",
)[
  Replace $V_(0 , h) := S_p^0 (cal(M)) , p in bb(N)$, with
  $V'_(0 , h) := S_(p + 1)^0 (cal(M))$, yielding a larger space: $V_(0 , h) subset V'_(0 , h)$

]
So h-refinement refines the mesh (smaller and smaller triangles) and
p-refinement chooses more powerful basis functions (start with linear, then
quadratic, etc.). The h in h-refinement refers to the mesh width:
#definition(
  number: "3.2.1.4", "Mesh width",
)[
  Given a mesh $cal(M) = { K }$, the *mesh width* $h_(cal(M))$ is defined as
  $ h_(cal(M))        & := max { upright("diam") K : K in cal(M) }\
  upright("diam") K & := max { lr(|p - q|) : p , q in K } $
]

== Empirical (Asymptotic) Convergence of Lagrangian FEM
<sub:empirical-asymptotic-convergence>

As in NumCSE, there are two types of convergence, algebraic and exponential. We
refer to the number of basis functions (dimension of the trial space) as $N$ and
we study the behavior of errors as
$N arrow.r oo$. Note that both exercises and theorems are often posed in terms
of $h_(cal(M))$ as $h arrow.r 0$, which is equivalent if we have a fixed
polynomial degree $p$ and only do h-refinement.

#definition(
  number: "3.2.2.1", "Types of convergence",
)[

  $ norm(u - u_N)_a = Order(N^(- alpha)) , alpha > 0$ is called *algebraic*
  convergence with rate $alpha$. \

  $ norm(u - u_N)_a = Order(exp (- gamma N^delta)) , gamma , delta > 0$ is called
  *exponential* convergence. ]

#tip-box(
  [Determining convergence rates], [
    - *algebraic*
    $ alpha approx (log eps_(i-1) - log eps_i)/(log N_i - log N_(i-1)) = log(eps_i\/eps_(i-1))/log(h_i\/h_(i-1)) $
    - *exponential* In general: complicated (see §3.2.2.5 Lecture Notes)

      If $delta=1$ (plain exponential convergence) and we let $N$ increase linearly
      (e.g., $N_i = i$), $ eps_(i+1)/eps_i approx exp(-gamma) $
  ],
)
Note that in the case of h-refinement we get the relation between $N$
and $h$ given by
#equation(
  number: "3.2.2.1",
)[
  #neq(
    $ N = dim S_p^0 (cal(M)) approx p^d h_(cal(M))^(- d) thick arrow.r.double.long h_(cal(M)) / p approx N^(- 1 / d) $,
  ) <eq:h-n-relation>
  where $p$ are the dimensions of the local basis functions, e.g., $p = 1$
  for linear basis functions, and $d$ is the dimension of the underlying space $Omega$.
]

E.g., in the case where we have $Omega subset bb(R)^2$ and piecewise linear
basis functions, we get $h_(cal(M))^(- 2) approx N$.

== A Priori (Asymptotic) Finite Element Error Estimates
<sub:a-priori-asymptotic-error>

Since FEM is similar to polynomial interpolation, we can use interpolation error
estimates to get error bounds. Here are some results for linear interpolation:
#mybox(
  "Linear interpolation error 1D",
)[
  Using the linear interpolant $I_1$ we want to study the interpolation error $u - I_1 u$.
  The following interpolation error estimates can be used for sufficiently smooth
  functions $u$:
  #neq(
    $ norm(u - I_1 u)_(L^oo (openint(a, b))) <= 1 / 4 h_(cal(M))^2 norm(u'')_(L^oo (openint(a, b))) $,
  )
  #neq(
    $ norm(u - I_1 u)_(L^2 (openint(a, b))) <= h_(cal(M))^2 norm(u'')_(L^2 (openint(a, b))) $,
  )
  #neq(
    $ lr(|u - I_1 u|)_(H^1 (openint(a, b))) <= h_(cal(M)) norm(u'')_(L^2 (openint(a, b))) $,
  )
]
#mybox(
  "Linear interpolation error 2D",
)[
  In 2D, linear interpolation corresponds to using tent functions: $I_1 u = sum_(p in cal(V) (cal(M))) u (p) b^p$,
  where $b^p$ is the tent function associated with point $p$.
  #neq(
    $ norm(u - I_1 u)_(L^2 (Omega))        &<= C h_(cal(M))^2 norm(norm(D^2 u)_F)_(L^2 (Omega)) \
    norm(grad (u - I_1 u))_(L^2 (Omega)) &<= C rho_(cal(M)) h_(cal(M)) norm(norm(D^2 u)_F)_(L^2 (Omega)) $,
  ) <eq:lin_interpolate_2d>
  Here, $D^2 u$ is the Hessian of $u$, $norm(dot)_F$ is the Frobenius norm, and $rho_(cal(M))$ is
  the shape regularity measure of the mesh $cal(M)$, defined as $rho_(cal(M)) = max_(K in cal(M)) h_(cal(M))^2 / lr(|K|)$ for
  a triangular mesh.
]

To get rid of this cumbersome notation, we can introduce more Sobolev spaces.

#definition(
  number: "3.3.3.1", "Higher order Sobolev spaces/norms", ..unimportant,
)[
  The $m$-th order Sobolev norm is defined as
  #neq(
    $ norm(u)_(H^m (Omega))^2 = sum_(k = 0)^m sum_(balpha in bb(N)^d , lr(|balpha|) = k) integral_Omega lr(|D^balpha u|)^2 dif bx, wide "where" D^balpha u = frac(
      partial^(lr(|balpha|)) u, partial x_1^(alpha_1) dots.h.c thin partial x_d^(alpha_d),

    ) $,
  )
  Hence we can define the $m$-th Sobolev space as
  #neq(
    $ H^m (Omega) = {v : Omega arrow.r bb(R) : norm(v)_(H^m (Omega)) < oo} $,
  )
]
#definition(
  number: "3.3.3.3", "Higher order Sobolev semi-norms", ..unimportant,
)[
  The $m$-th order Sobolev semi-norm is defined as
  #neq(
    $ lr(|u|)_(H^m (Omega))^2 = sum_(bold(alpha) in bb(N)^d , lr(|bold(alpha)|) = m) integral_Omega lr(|D^bold(alpha) u|)^2 dif bx $,
  )
]

Remember the multi-index $balpha$ already seen in @concept:multi-index. Using
this new notation, we can rewrite the error bounds from @eq:lin_interpolate_2d
as
$ norm(u - I_1 u)_(L^2 (Omega))        &<= C h_(cal(M))^2 lr(|u|)_(H^2 (Omega)) \
norm(grad (u - I_1 u))_(L^2 (Omega)) &<= C rho_(cal(M)) h_(cal(M)) lr(|u|)_(H^2 (Omega)) $
It turns out that these bounds are not sharp. There is a very useful result for
Lagrangian finite elements:
#theorem(
  number: "3.3.5.6", "Best approximation error estimates for Lagrangian finite elements",
)[
  Given a triangular mesh $cal(M)$, if the true solution $u$ is in $H^k (Omega)$,
  the best approximation error is bounded by
  #neq(
    $ inf_(v_h in cal(S)_p^0 (cal(M))) norm(u - v_h)_(H^1 (Omega)) <= C (h_(cal(M)) / p)^(min { p , k - 1 }) norm(u)_(H^k (Omega)) $,
  )
] <thm:best-approximation-error>

We might not know the constant $C$ and/or $norm(u)_(H^k (Omega))$, but we know $p$ and $h_(cal(M))$ as
they are imposed by the choice of function space and mesh. Remember the concept
of refinement: we can adjust these values. And from Eq. @eq:h-n-relation we know
that $h_(cal(M)) \/ p approx N^(- 1 / d)$. Hence the error displays *algebraic*
convergence with rate
$min{ p , k - 1 }\/ d$. What still remains a question is $k$, the smoothness of
the solution $u$.

Another useful result hidden in the lecture notes is
#equation(
  number: "3.6.3.10", [$L^2$ estimate],
)[
  Under some assumptions (convex domain, smooth coefficient functions), we have
  #neq(
    $ norm(u - u_h)_(L^2 (Omega)) <= C h_(cal(M)) / p norm(u - u_h)_(H^1 (Omega)) $,
  )
]
So we gain one order of convergence in the $L^2$ norm compared to the $H^1$ norm.
#tip-box(
  [Rules of thumb for converge], [
    If we are using $cal(S)_p^0 (cal(M))$ and $u$ is sufficiently smooth (e.g., $u in C^oo (Omega)$),
    we have
    $ norm(u - u_h)_(H^1 (Omega)) &= Order(h^p) \
    |u - u_h|_(H^1 (Omega))     &= Order(h^p) \
    norm(u - u_h)_(L^2 (Omega)) &= Order(h^(p + 1)) $
  ],
)

== Elliptic regularity
<sub:elliptic-regularity>

#theorem(
  number: "3.4.0.2", "Smooth elliptic lifting theorem",
)[
  For domains $Omega$ with smooth boundaries $partial Omega$, i.e. no corners and
  sufficiently smooth $sigma$, if
  $ u in H_0^1 (Omega) quad upright("and") quad - div (sigma grad u) in H^k (Omega) $
  or
  $ u in H^1 (Omega) , - div (sigma grad u) in H^k (Omega) quad upright("and") quad grad u dot.op bold(n) = 0 quad upright("on") thin partial Omega $
  holds, then $u in H^(k + 2) (Omega)$ and
  $ norm(u)_(H^(k + 2) (Omega)) <= C norm(div (sigma grad u))_(H^k (Omega)) $
]

This tells us that when solving the typical PDE
$-div (sigma grad u) = f$ and the source term $f$
is in $H^k (Omega)$, the solution $u$ will be in $H^(k + 2) (Omega)$ (of course
under the right assumptions).

The theorem requires smooth domains, but our meshes will have corners, so what
can be done there? As long as the domain and all cells are convex, something
similar still holds:
#theorem(
  number: "3.4.0.10", "Elliptic lifting on convex domains",
)[
  If $Omega subset RR^d$ is convex, $u in H_0^1 (Omega)$ and $Delta u in L^2 (Omega)$,
  then $u in H^2 (Omega)$.
]<thm:elliptic-lifting-convex-domains>
If we are solving the Laplace equation on a convex domain, we just need to check
if $f in L^2 (Omega)$, since $-Delta u = f$.
#tip-box(
  [Finding $k$ in @thm:best-approximation-error], [
    - Often, the exercise gives you $u$ which is in $C^oo (Omega)$, which means $u$ is
      infinitely smooth. In this case, $k = oo$.
    - Sometimes, you may use elliptic regularity results like
      @thm:elliptic-lifting-convex-domains to find $k=2$, for example.
  ],
)

== Variational Crimes
<sub:variational-crimes>
What are variational crimes?

A variational crime is committed when instead of the true variational problem
$ u_h in V_(0 , h) : thin a (u_h , v_h) = ell (v_h) , thin forall v_h in V_(0 , h), $
we solve a different, "perturbed" variational problem
$ tilde(u)_h in V_(0 , h) : thin a_h (tilde(u)_h , v_h) = ell_h (v_h) , thin forall v_h in V_(0 , h) $
with modified (bi-)linear forms $a_h , ell_h$. With computers, the use of
quadrature and approximation of boundaries result in such a crime and are
unavoidable. As Hiptmair likes to say, "we are all sinners".

So the only distinction we can make is between acceptable and unacceptable "crimes".
Crimes which do not affect the type and rate of convergence are acceptable.

So how to not temper with the convergence?

#subtle-box[
  - if $norm(u - u_h)_1 in Order(h_(cal(M))^p)$, use a quadrature rule of order at
    least $2 p - 1$

  - if $V_(0 , h) = cal(S)_p^0 (cal(M))$ then approximate the boundary with
    polynomials of degree $p$
]

#pagebreak()
== FEM: Duality Techniques for Error Estimation
<sub:duality>

#theorem(
  number: "3.6.1.7", "Duality estimate for linear functional output",
)[
  Given a functional $F : V_0 arrow.r bb(R)$ the dual solution $g_F$ solves
  $ g_F in V_0 : thin a (g_F , v) = F (v) , thin forall v in V_0 $
  and we get the estimate
  $ lr(|F (u) - F (u_h)|) <= C norm(u - u_h)_a inf_(v_h in V_(0 , h)) norm(g_F - v_h)_a $
]

Why is this useful? If $g_F$ can be approximated well in $V_(0 , h)$, then the
output error $lr(|F (u) - F (u_h)|)$ can converge to 0 much faster than $norm(u - u_h)_a$.