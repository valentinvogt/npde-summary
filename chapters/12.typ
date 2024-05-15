#import "../src/setup.typ": *
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
$ bv^* &= argmin_(bv in V) &1/2 integral_Omega mu norm(curl bv)^2 dif bx - integral_Omega bold(f) dot bv dif bx \
&= argmin_(bv in V) &1/2 integral_Omega mu norm(jac bv)_F^2 dif bx - integral_Omega bold(f) dot bv dif bx
$
where $bold(f)$ is a given force field, $mu$ is the _viscosity_.

We also have an equivalent variational (weak) form
$ bv^* in V upright("such that ") wide integral_Omega mu med jac bv : jac bw dif bx = integral_Omega bold(f) dot bw dif bx wide forall bw in V $