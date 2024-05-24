#import "../setup.typ": *
#show: thmrules

= Finite Elements for the Stokes Equation
<ch:stokes>

#counter(heading).update((12, 2))
=== Constrained Variational Formulation
To simulate a fluid, we solve for a velocity field $bv$, so we have a vector
problem. We assume incompressibility (recall @sub:heat-conduction) and the _non-slip boundary condition_ $bv = 0$ on $partial Omega$.
The Stokes equation is then
$ V = {bv : Omega -> RR^d upright("continuous"), div bv = 0, bv = 0 upright("on") partial Omega} $

We basically generalize @ch:01 to vector fields:
$ u                    &--> bv \
grad u               &--> jac bv          &&upright("(Jacobian)") \
norm(grad u)         &--> norm(jac bv)_F  &&upright("(Frobenius norm)") \
grad u dot.op grad v &--> jac bv : jac bw && \ $
where the component-wise dot product is defined as $A : B= sum_(i j) A_(i j) B_(i j)$.
We also get a Sobolev space
$ bH^1_0(div 0,Omega) = {bv in (H^1_0(Omega))^d : div bv = 0} $

Like in @ch:01, we have a minimization problem
#neq(
  $ bv^* &= argmin_(bv in V) &1/2 integral_Omega mu norm(curl bv)^2 dif bx - integral_Omega bold(f) dot bv dif bx \
       &= argmin_(bv in V) &1/2 integral_Omega mu norm(jac bv)_F^2 dif bx - integral_Omega bold(f) dot bv dif bx $,
)<eq:stokes-minimization>
where $bold(f)$ is a given force field, $mu$ is the _viscosity_.

There is a variational (weak) form equivalent to @eq:stokes-minimization:
#neq(
  $ bv^* in V upright("such that ") wide integral_Omega mu med jac bv : jac bw dif bx = integral_Omega bold(f) dot bw dif bx wide forall bw in V $,
)<eq:stokes-variational>
It's difficult to solve this directly because of the constraint $div bv = 0$.
Therefore, we use the method of Lagrange multipliers.

#pagebreak()

=== Saddle Point Formulation

We apply the method of Lagrange multipliers (generalized to vector fields) to
@eq:stokes-variational to get the following "saddle point" problem:

#subtle-box(
  )[
  Seek the velocity field $bv in (H_0^1(Omega))^d$ and a Lagrange multiplier $p in L^2_*(Omega)$ such
  that
  #neq(
    $   &integral_Omega mu jac bv : jac bw dif bx + integral_Omega div bw med p dif bx &&= integral bold(f) dot bw dif bx quad &&forall bw &&in (H_0^1(Omega))^d\
      &integral_Omega div bv med q dif bx                                            &&= 0 quad                              &&forall q  &&in L^2_*(Omega) $,
  )<eq:ssp>
]
We can interpret the Lagrange multiplier $p$, which is a scalar field, as the _pressure_.

As trial and test space for the pressure, we use the space
#neq($ L^2_* = {q in L^2(Omega) : integral_Omega q dif bx = 0} $)<eq:l2-star>
because like in @sub:second-order-elliptic-variational-problems, the pressure is
only defined up to a constant (if $p$ solves @eq:ssp, then $p + c$ also solves
it for any $c in RR$).

In practice, to deal with this constraint, we enforce it through a Lagrange multiplier too:
#subtle-box[
  Seek $bv in (H_0^1(Omega))^d$, $p in L^2(Omega)$, and $lambda in RR$ such that
  #neq(
    $
    & integral_Omega mu jac bv : jac bw dif bx + && integral_Omega div bw med p dif bx & quad & quad =  integral bold(f) &&dot bw dif bx quad &&forall bw in (H_0^1(Omega))^d \
    & integral_Omega div bv med q dif bx && quad & + lambda integral_Omega q dif bx & quad = && 0 quad &&forall q in L^2(Omega) \
    & quad && integral_Omega p dif bx & &quad = && 0
    $
)<eq:ssp-zero-mean>
]

As always, we are concerned with the existence and uniqueness of solutions to
@eq:ssp.
#theorem(
  number: "12.2.2.40", [Existence and uniqueness of weak solutions of Stokes problem],
)[
  The Stokes problem @eq:ssp has a unique solution $(bv, p) in bold(H)^1_0(div 0,Omega) times L^2_*(Omega)$ which
  satisfies
  #neq(
    $ norm(bv)_(H^1(Omega)) + norm(p)_(L^2(Omega)) <= C norm(bold(f))_(L^2(Omega)) $,
  )<eq:stokes-estimate>
  where $C= C(Omega)$.
]

#pagebreak()
=== Stokes System of PDEs

The strong form of @eq:ssp is the following system of PDEs:
#align(center)[
  #subtle-box(width: 50%)[
    #v(-0.2cm)
    $ -mu bold(Delta) bv + grad p &= bold(f) \
    div bv                      &= 0 upright("on ") Omega \
    integral_Omega p dif bx     &= 0 \
    bv                          &= 0 upright("on ") partial Omega $
    #v(-0.2cm)
  ]
]

Here, $bold(Delta)$ is the component-wise Laplacian.

== Galerkin Discretization of the Stokes Equation
<ch:stokes-galerkin>

We can write the saddle point problem @eq:ssp in a more abstract form: 

Let $U= (H^1_0(Omega))^d$ and $Q= L^2_*(Omega)$.
Seek $bv in U$ and $p in Q$ such that
#neq($ a(bv, bw) + &b(bw, p) &&= ell(bw) && quad forall bw in U \
                   &b(bv, q) &&= 0       && quad forall q in Q $)
                   <eq:ssp-abstract>
where
$ a(bv, bw) = integral_Omega mu jac bv : jac bw dif bx , wide
b(bv, q) = integral_Omega div bv med q dif bx $

To apply the Galerkin method, we need basis functions for discrete subspaces $U_h subset U$ and $Q_h subset Q$:
$ frak(B)_U = {bold(phi)_h^1, ..., bold(phi)_h^N}, wide
  frak(B)_Q = {beta_h^1, ..., beta_h^M} $

where $N= dim U_h$ and $M= dim Q_h$. The basis of $U_h$ is now vector-valued, but we can simply use scalar-valued basis functions for each component, e.g. $U_h = (cal(S)_1^0)^d$.

With a basis expansion, we can define Galerkin matrices:
$ bv_h = sum_(i=1)^N nu_i bold(phi)_h^i wide bold(arrow(nu))=[nu_i]_(i=1)^N \
  p_h  = sum_(i=1)^M pi_i beta_h^i wide bold(arrow(pi))=[pi_i]_(i=1)^M $
$ bold(A) = [a(bold(phi)_h^j, bold(phi)_h^i)]_(i, j=1)^N wide
  bold(B) = [b(bold(phi)_h^j, beta_h^i)]_(i, j=1)^(M, N) wide
  bold(arrow(gamma)) = [ell(bold(phi)_h^i)]_(i=1)^N $

And we get the following linear system:
$ upright("Saddle point LSE: ") quad 
bnat(
  bold(A), bold(B)^T;
  bold(B), bold(0)
) bnat(
  bold(arrow(nu)); bold(arrow(pi))
) = bnat(
  bold(arrow(gamma)); bold(0)
) $


=== Pressure instability
<sub:stokes-pressure-instability>

It turns out that not every choice of discrete spaces leads to a stable method.

For the pressure field, we also consider _discontinuous_ piecewise polynomials from $L^2_*(Omega)$, see @eq:l2-star:
$ S^(-1)_(p,*)(msh) = {q_h in L^2(Omega) : quad eval(q_h)_K in cal(P)_p (RR^d) , quad integral_Omega q_h dif bx = 0 } $

Now we can analyze methods with respect to stability and convergence.

To study the convergence of the Galerkin method, we consider the error 
$ E(bv_h, p_h) := norm(bv - bv_h)_(H^1(Omega)) + norm(p - p_h)_(L^2(Omega)) $

The following theorem states that stable methods are _quasi-optimal_ (cf. @thm:ceas-lemma, where the Galerkin solution was optimal).
#theorem(
  number: "12.3.3.13", [Convergence of stable FE for Stokes problem],
)[
 If the method defined by $U_h$, $Q_h$ is stable, then the error satisfies
  $ E(bv_h, p_h) <= C (inf_(bw_h in U_h) norm(bv - bw_h)_(H^1(Omega)) + inf_(q_h in Q_h) norm(p - q_h)_(L^2(Omega))) $
  where $C$ is a constant depending only on the mesh and $mu$.
]

We can apply our error estimates from @ch:convergence-and-accuracy to this to get convergence rates. For example, take $U_h = (cal(S)_1^0)^d$ and $Q_h = S^(-1)_(0,*)(msh)$, the so-called *P1-P0* method. 

Then we have
$ E(bv_h, p_h) <= C h^2 norm(bv)_(H^3(Omega)) + C h norm(p)_(H^1(Omega)) = Order(h^2) + Order(h) = Order(h) 
$ The first terms comes from @thm:best-approximation-error assuming $bv$ to be smooth, the second term comes from @eq:l2-estimate.

Some methods for the Stokes problem are:
#table(
  columns: 6,
  inset: 7pt,
  align: (left, center, center, left, center, center),
  [Method], [$U_h$], [$Q_h$], [Stability], [Convergence], [efficient?],
  [*P1-P0*], [$(cal(S)_1^0(msh))^d$], [$S^(-1)_(0,*)(msh)$], [Unstable], [$Order(h)$], [$checkmark$],
  [*P2-P0*], [$(cal(S)_2^0(msh))^d$], [$S^(-1)_(0,*)(msh)$], [Stable], [$Order(h)$], [$times$],
  [*P2-P1*], [$(cal(S)_2^0(msh))^d$], [$S^0_(1,*)(msh)$], [Stable], [$Order(h^2)$], [$checkmark$],
  [*CR-P0*], [$(cal(C R)_0(msh))^d$], [$S^(-1)_(0,*)(msh)$], [Stable], [$Order(h)$], [$checkmark$], 
)

The $cal(C R)$ space will be introduced in the next subsection. By "efficient" we mean that both error terms convergence with the same rate. If that is not the case -- as for the *P2-P0* method, where the pressure converges slower -- the computational cost for having 2nd degree polynomials for the velocity does not pay off.

=== Non-conforming Crouzeix--Raviart FEM
<sub:stokes-cr>

We define a new discrete space of basis functions for the velocity field: $cal(C R)(msh)$, the _Crouzeix--Raviart_ space. It has one *discontinuous, piecewise linear* basis function associated with every *edge* of the mesh. Let $bm_e$ be the midpoint of edge $e$. The basis functions are defined as

$ eval(b^e_h)_K in cal(P)_1(K) &quad forall K in msh\ quad b^e_h (bm_f) = delta_(e f) &quad forall e,f in cal(E)(msh) $

Note that, in contrast to the Lagrangian spaces $cal(S)_p^0(msh)$, the Crouzeix--Raviart space is not a subspace of $H^1(Omega)$. This is why it is called _non-conforming_.