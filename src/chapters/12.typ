#import "../setup.typ": *
#show: thmrules

= Finite Elements for the Stokes Equation
<ch:stokes>
To simulate a fluid, we solve for a velocity field $bv$, so we have a vector
problem. We assume incompressibility (recall @sub:heat-conduction) and the _non-slip boundary condition_ $bv = 0$ on $partial Omega$.
The Stokes equation is then
$ V = {bv : Omega -> RR^d upright("continuous"), div bv = 0, bv = 0 upright("on") partial Omega} $

We basically generalize @ch:01 to vector fields:
$ u                    &--> bv \
grad u               &--> jac bv         &&upright("(Jacobian)") \
norm(grad u)         &--> norm(jac bv)_F &&upright("(Frobenius norm)") \
grad u dot.op grad v &--> jac bv : jac bw  && \ $
where the component-wise dot product is defined as $A : B= sum_(i j) A_(i j) B_(i j)$.
We also get a Sobolev space
$ bH^1_0(div 0,Omega) = {bv in (H^1_0(Omega))^d : div bv = 0} $

Like in @ch:01, we have a minimization problem
#neq($ bv^* &= argmin_(bv in V) &1/2 integral_Omega mu norm(curl bv)^2 dif bx - integral_Omega bold(f) dot bv dif bx \
&= argmin_(bv in V) &1/2 integral_Omega mu norm(jac bv)_F^2 dif bx - integral_Omega bold(f) dot bv dif bx
$)<eq:stokes-minimization>
where $bold(f)$ is a given force field, $mu$ is the _viscosity_.

There is a variational (weak) form equivalent to @eq:stokes-minimization:
$ bv^* in V upright("such that ") wide integral_Omega mu med jac bv : jac bw dif bx = integral_Omega bold(f) dot bw dif bx wide forall bw in V $
Now it is hard to solve this directly because of the constraint $div bv = 0$. Therefore, we use the method of Lagrange multipliers.

#counter(heading).step(level: 2)
== Saddle Point Formulation
// Seek the velocity field $\mathbf{v}\in(H_0^1(\Omega))^d$ and a Lagrange multiplier $p\in L^2(\Omega)$ such that
// $$\begin{aligned}&\int\mu\mathrm{D}\mathbf{v}:\mathrm{D}\mathbf{w}\:\mathrm{d}\mathbf{x}\:+\:\int\mathrm{d}\mathrm{i}\mathrm{v}\:\mathbf{w}\:p\:\mathrm{d}\mathbf{x}\:=\:\int\mathbf{f}\cdot\mathbf{w}\:\mathrm{d}\mathbf{x}\:\quad\forall\mathbf{w}\in(H_0^1(\Omega))^d\:,\\&\int_\Omega\operatorname{div}\mathbf{v}\:q\:\mathrm{d}\mathbf{x}\:=\:0&\forall q\in L^2(\Omega)\:.\end{aligned}$$

#subtle-box()[
  Seek the velocity field $bv in (H_0^1(Omega))^d$ and a Lagrange multiplier $p in L^2(Omega)$ such that
  $ &integral_Omega mu jac bv : jac bw dif bx + integral_Omega div bw med p dif bx &&= integral bold(f) dot bw dif bx quad &&forall bw &&in (H_0^1(Omega))^d\
  &integral_Omega div bv med q dif bx &&= 0 quad &&forall q &&in L^2(Omega) $
]

This problem has the following strong form:
#subtle-box()[
  $ -mu Delta bv + grad p &= bold(f) \
  div bv &= 0 upright("on ") Omega \
  integral_Omega p dif bx &= 0 \
  bv &= 0 upright("on ") partial Omega $
]
Here, $Delta$ is the component-wise Laplacian.
