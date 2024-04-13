#import "../src/setup.typ": *
#show: thmrules

= Non-Linear Elliptic Boundary Value Problems
<ch:non-linear-elliptic-bvp>
== Elastic String Model
<sub:elastic-string-model>
We want to derive the general variational equation for an elastic
string. For this, one approximates the string as $n$ point masses
affected by gravity connected with springs, whose energy behaves
according to Hooke’s law. Then, one takes the limit $n arrow.r oo$ to
derive a continuous model. The total energy is then just given by the
sum of elastic and gravitational energies — given positions
$(mu_0 , dots.h , mu_n)$ and $x_i = a + h i$, assuming the spring
constants are $1$:
#mybox("Total energy of discrete spring system", ..unimportant)[
  $ E (mu) = 1 / 2 sum_(i = 0)^n (sqrt(h^2 + (mu_(i + 1) - mu_i)^2))^2 + sum_(i = 1)^n m_i mu_i g $
]
Then, the equilibrium position for this model can be found by minimizing
this expression over $mu$. A continuous model is derived by replacing
the discrete positions $mu_i$ by a function $u (x_i)$, and the mass by a
mass density. Then, performing some manipulations, one obtains

#mybox("Total energy for the continuous string model")[
  $ J_s (u) = integral_a^b 1 / 2 frac(b - a, L) sigma (x) (sqrt(1 + lr(|u prime (x)|)^2) - frac(L, b - a))^2 + integral_a^b g rho (x) u (x) dif x $
]
where $sigma$ is the spring stiffness, $L$ is the total rest length of the string, and $rho$ is the mass density. The term $sqrt(1 + lr(|u' (x)|)^2)$ is the length of the spring at position $x$.

In a similar fashion, a membrane model can be derived by assuming a
two-dimensional grid of springs containing points masses and taking the
limit $n arrow.r oo$ springs. The energy then becomes
#mybox("Total energy for the membrane model", ..unimportant)[
  $ J_M (u) = integral_Omega frac(1, 2 L) sigma (bx) (&(sqrt(1 + lr(|frac(partial u, partial x_1) (bx)|)^2) - frac(L, b - a))^2 \
  + &(sqrt(1 + lr(|frac(partial u, partial x_2)(bx)|)^2) - frac(L, b - a))^2) + g rho (bx) u (bx) dif bx $

]
In the limit of a taut membrane, i.e. for $L lt.double b - a$, these
equations just reduce to the problem of minimizing the familiar
functionals seen in earlier chapters:

#mybox("Elastic string/membrane: taut limit")[
  #neq($ J_s (u) &= 1 / 2 integral_a^b hat(sigma) (x) lr(|u prime (x)|)^2 + g rho (x) u (x) dif x \
  J_M (u) &= 1 / 2 integral_Omega hat(sigma) (bx) norm(grad med u (bx))^2 + g rho (bx) u (bx) dif bx $)<eq:taut-limit>
]

== Calculus of Variations
<sub:calculus-of-variations>
The difference between the equations seen in @eq:taut-limit, which hold in the
limit of a very stretched string/membrane, and the general equations
given above, is that the former are #strong[quadratic] minimization
problems, while the latter are #strong[nonlinear] minimization problems.
The theory used so far mapped quadratic minimization problems to
#strong[linear] variational problems, which were then discretized. The
new equations, however, yield nonlinear variational equations.
Therefore, more general variational problems need to be derived.

The idea employed is that, for a minimizer $u_*$ of $J (u)$, every
perturbation $J (u_* + v)$ would be larger than $J (u_*)$. This means that
$f (t) = J (u_* + t v)$ has a minimum at $t = 0$ for every function $v$:

#theorem(number: "5.2.1.5", "Characterization of global minimizers")[
  Assume $u_* in V_0$ is a global minimizer of $J (u)$, i.e. $ u_* = op("argmin", limits: #true)_(u in V_0) J (u) $ Then, if $phi_v (t) = J (u_* + t v)$ is differentiable in $t = 0$, we have $ frac(dif phi_v, phi t) (0) = 0 quad forall v in V_0 $
]

This means that nonlinear variational equations can be derived by
computing this derivative for an arbitrary $v$. As an example, for the
elastic string model introduced in the last sub-chapter, this yields
#mybox("Variational equations for elastic string model")[
  $ integral_a^b frac(sigma (x), c) (sqrt(1 + lr(|u prime (x)|)^2) - c) frac(u prime (x) v (x), sqrt(1 + lr(|u prime (x)|)^2)) + g rho (x) v (x) dif x = 0 quad forall v in H_0^1 (\] a , b \[) $
]
This can be formulated more generally as a #strong[general variational equation];:
#mybox("General variational equation")[
  A general, nonlinear variational equation reads
  $ u in mhat(V) : quad a (u ; v) = 0 quad forall v in V_0 $
  Where $a$ is #strong[linear in the second argument] $v$ and
  $V_0 , mhat(V)$ are function spaces.
]

== Nonlinear Boundary Value Problems
<sub:nonlinear-bvp>
Similarly to the linear case, nonlinear PDEs can be derived from the
variational equations by "stripping" of the derivatives of $v$ by
partial integration and employing the fundamental lemma of calculus of
variations. As an example, for the string model, this yields for
$u (a) = u_a$, $u (b) = u_b$

$ frac(dif, dif x) (frac(sigma (x), c) (sqrt(1 + lr(|u '(x)|)^2) - c) frac(u' (x), sqrt(1 + lr(|u ' (x)|)^2))) = g rho (x) quad upright(" in ") thin \] a , b \[ $

== Galerkin Discretization of Non-Linear BVPs
<sub:galerkin-non-linear-bvps>
The idea of Galerkin discretization for non-linear variational equations
is exactly the same as for linear equations, but they yield nonlinear
systems of equations instead of linear systems of equations: One
restricts $u$ and $v$ to a finite function space $u_h in mhat(V)_h$ and
$v_h in V_(0 , h)$, and expands the functions in some basis of the
space. This, then, leads to nonlinear equations for the basis expansion
coefficients. These equations could be solved directly by employing some
fixed-point iteration seen in NumCSE.
#mybox("Nonlinear Galerkin Discretization")[
  Given a variational problem $a (u ; v) = 0 med forall v in V_(0 , h)$, the Galerkin discretization reads
  $ (F (mu))_i = a (u_(0 , h) + sum_(j = 1)^N mu_j b_h^j ; b_h^i) , quad i = 1 , dots.h , N $
  where $b_h^i$ are fixed basis functions, $mu_i$ are the basis function coefficients and $u_(0 , h)$ contains Dirichlet boundary conditions.
]
Another option is to already linearize the continuous problem, and then
discretize it to derive linear systems of equations. This is done by
employing Newton's method in function space: The conventional Newton
iteration is given as
$ xi^((k + 1)) = xi^((k)) - upright(D) F (xi^((k)))^(- 1) F (xi^((k))) $
Replacing the vector $xi$ with a function $u$ and the derivative by a
function derivative now gives
#equation(number: "5.3.2.5", "Functional Newton iteration")[
  $ w in V_0 : a (u^((k)) ; v) + upright(D)_u a (u^((k)) ; v) w = 0 quad forall v in V_0\
  u^((k + 1)) = u^((k)) + w $
]

Here, the directional derivative is defined as
#align(center)[
  #subtle-box[
    $ upright(D)_u a (u^((k)) ; v) w = lim_(t arrow.r 0) frac(a (u + t w ; v) - a (u ; v), t) , quad u^((k)) in mhat(V) , quad v , w in V_0 $
  ]
]
#tip("Computing the directional derivative")[
  You can simply compute the defining limit by Taylor-expanding $a(u + t w, v)$. If you do this once, you notice that only the first-order term in $t$ remains: the 0th order term cancels out with $a(u, v)$, and the higher-order terms vanish in the limit $t arrow.r 0$.
  
  So a more convenient method is to do Taylor expansion and to keep only the first-order term.
]
Now, the advantage of this equation for $w$ is that the functional
derivative is linear, i.e. $(v , w) arrow.r.bar D_u a (u^((k)) ; v) w$
is a #strong[bilinear form];. Now, one can employ Galerkin
discretization for the linear problem in $w$, exactly like it was done
in Chapter 2 and 3. The final equations then read

#equation(number: "5.3.3.6", "Nonlinear Newton equations for variational problems")[
  $ w_h in V_(0 , h)^((k)) : upright(D)_u a (u_h^((k - 1)) ; v_h) w_h = - a (u_h^((k - 1)) ; v_h) quad forall v_h in V_(0 , h)^((k))\
  u_h^((k)) = P_h^((k)) (u_h^((k - 1)) + w_h) $
]
Here, different function spaces can be used for each iterations, so a
projector $P_h^((k))$ needs to be used to project the solution from
$V_h^((k - 1))$ to $V_h^((k))$. In all of these equations, the previous
iterate $u_h^((k - 1))$ is kept fixed, a linear system like derived in
Chapter 2 is solved to obtain the intermediate $w_h$, and then a new
iterate $u^((k))$ is obtained.
#pagebreak()