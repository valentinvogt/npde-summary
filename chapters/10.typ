#import "../src/setup.typ": *

#show: thmrules

= Convection-Diffusion Problems
<ch:convection-diffusion>
== Heat conduction in a Fluid
<sub:heat-conduction>
Consider a flowing fluid. Then there is the key quantity, the #emph[flow
field] $bold(v) : Omega subset bb(R)^d arrow.r bb(R)^d$, where $d$ is
the dimension we consider. The flow field can be understood as
$bold(v) (bx) =$ fluid velocity at point $bx in Omega$.

Given a flow field $bold(v)$, we can consider the autonomous initial
value problem
$ dot(bold(y)) = bold(v) (bold(y)) , quad bold(y)_0 = bx_0 $ The
solution $t arrow.r bold(y) (t)$ describes how a particle moves, carried
by the fluid, also called #emph[streamline];.

As the domain $Omega$ is usually bounded, we cannot have fluid leaving
the domain. This means the fluid velocity must be zero in the normal
direction of the domain boundary. Hence
$ bold(v (x) dot.op n (x)) = 0 , quad forall bx in partial Omega $


#strong[Fourier’s law in a moving fluid]
$ bold(j (x)) = - kappa grad  u (bx) + bold(v (x)) rho u (bx) $
with $kappa$ the heat conductivity and $rho$ the volumetric heat
capacity. We already know the first part, called diffusive heat flux, from the
heat equation. The second part is called convective heat flux. With this
new flux, the standard PDE becomes
#neq($ - div (kappa grad u) + div (rho bold(v (x)) u) = f quad upright("in ") Omega $) <eq:convection-diffusion-strong>
#strong[Incompressible Fluids]

A fluid is called incompressible if its
associated flow map (evaluation operator) $Phi^t$ is volume preserving. This is the case iff.
#align(center)[
  #subtle-box(width: 50%)[
    $ upright("div ") bold(v) equiv 0 quad upright("in ") Omega $ 
  ]
]
Hence the fluid is incompressible if its flow velocity is divergence-free. 

// #theorem(number: "10.1.3.7", "Differentiation formula for determinants")[
//   Let $bold(S)$ be a smooth matrix-valued function. If $bold(S)$ is
//   regular then
//   $ frac(dif, dif t) upright("det") bold(S) (t) = upright("det") (bold(S) (t)) upright("tr") (frac(dif bold(S), dif t) (t) bold(S)^(- 1) (t)) $
//   with det the determinant and tr the trace of a matrix.
// ]

In case of incompressibility, equation @eq:convection-diffusion-strong can be
simplified using the general product rule @thm:general-product-rule and
$upright("div ") bold(v) = 0$
$ - div (kappa grad u) + bold(v (x)) dot.op grad rho u = f quad upright("in ") Omega $
We can also look at the time-dependent heat flow, which is similar to
the heat equation @eq:heat-local
#neq($ frac(partial, partial t) med (rho u) - div (kappa grad u) + div (rho med bold(v) (bx ,t) med bu) = f quad upright("in ") Omega $)<eq:time-dependent-heat>
We will see later on how to solve this.

#pagebreak()
== Stationary Convection-Diffusion Problems
<sub:stationary-convection-diffusion>
Here we will focus on the convection diffusion equation
@eq:convection-diffusion-strong with constant $kappa$, $rho$, incompressible flow
$bold(v)$ and zero Dirichlet boundary conditions. Hence
$ - kappa Delta u + rho bold(v (x)) dot.op grad u = f quad upright("in ") Omega , #h(2em) u = 0 quad upright("on ") partial Omega $
Non-dimensionalizing the problem results in
$ - eps Delta u + bold(v (x)) dot.op grad u = f quad upright("in ") Omega , #h(2em) u = 0 quad upright("on ") partial Omega $
with $norm(bold(v))_(L^oo (Omega)) = 1$. This results in the following
variational form
#neq($ eps integral_Omega grad u dot.op grad w dif bx + integral_Omega (bold(v) dot.op grad u) w dif bx = integral_Omega f (bx) w dif x $) <eq:convection-diffusion-weak>
with the left-hand side the bilinear form $a (u , w)$. However, $a$ is
not symmetric. This also means that it does not induce an energy norm.
However, it is still positive definite (see lecture document).

=== Singular perturbation

A boundary value problem depending on a
parameter $eps$ is called singularly perturbed, if the limit
problem for $eps arrow.r eps_0$ is not compatible with
the boundary conditions.

For $eps = 0$, the above PDE is singular perturbed. It cannot
satisfy Dirichlet boundary conditions on the _outflow_ part of the
boundary. 
// TODO: give a clear example
$ Gamma_(upright("out")) &= bx in partial Omega : bold(v (x)) dot.op bold(n (x)) > 0 \
Gamma_(upright("in")) &= bx in partial Omega : bold(v (x)) dot.op bold(n (x)) < 0 $

=== Upwinding

When trying to solve Eq. @eq:convection-diffusion-weak with
the Galerkin approach, when $eps$ is very close to 0, one can
observe huge oscillations in the solution, which is not correct. It
comes from the fact that the Galerkin matrix becomes close to singular.
So our goal is to get a robust method that can solve Eq.
@eq:convection-diffusion-weak no matter the $eps$.

Consider again Eq. @eq:convection-diffusion-weak but in $d = 1$ and with zero boundary conditions
$ eps integral_0^1 frac(partial u, partial x) frac(partial w, partial x) dif x + integral_0^1 frac(partial u, partial x) w dif x = integral_0^1 f (x) w dif x $
To calculate the Galerkin matrix for an equidistant mesh with $M$ cells,
we use the global composite trapezoidal rule for the convective term
$ integral_0^1 psi (x) dif x approx h sum_(j = 0)^M psi (j h) $
Hence the convective term of the bilinear form will be approximated by
$ integral_0^1 frac(partial u_h, partial x) med w_h dif x approx h sum_(j = 1)^(M - 1) colMath(frac(partial u_h, partial x)  (j h), accentcolor) med w_h (j h) $
But $frac(partial u_h, partial x)  (j h)$ is not defined: $u_h in cal(S)_(1 , 0)^0$ is piecewise continuous, so its derivative has jumps at the nodes. However,
convection transports the information in the direction of $bold(v)$ ($= 1$
in our case). So we use the value of $frac(partial u_h, partial x)$ from "upwind", i.e., in direction $-bold(v)$, which here is the value in the previous cell:
$ frac(partial u_h, partial x)  (j h) = lim_(delta arrow.r 0) frac(partial u_h, partial x)  (j h - delta bold(v)) = eval(frac(partial u_h, partial x)) _openint(x_(j - 1) , x_j) $ 
And generalized in more dimensions
$ bold(v (p)) dot.op grad u_h (bold(p)) = lim_(delta arrow.r 0) bold(v (p)) dot.op grad u_h (bold(p) - delta bold(v (p))) $

#counter(heading).update((10,2,2,1))
==== Streamline Diffusion

A totally different idea for fixing the
problem of $eps arrow.r 0$ is to add some $h$-dependent
diffusion. I.e., replace $eps arrow.l eps + c (h)$ with
$c (h) > 0$. However, there is smearing in the internal layers. // TODO: Unclear
But the solution is smooth along the direction of $bold(v)$, so adding
diffusion along the velocity should not do any harm.

The method of #emph[Anisotropic diffusion] is born. On cell $K$, replace
$eps arrow.l eps bold(I) + delta_K bold(v)_K bold(v)_K^top$.
Here, $bold(v)_K$ is the local velocity, obtained by averaging over the vertices, and
$delta_K > 0$ some controlling parameter.

$ integral_Omega eps bold(I) + delta_K bold(v)_K bold(v)_K^top) grad u dot.op grad w dif bx + integral_Omega (bold(v) dot.op grad u) w dif bx = integral_Omega f (bx) w dif bx $
However this affects the solution $u$, such that it will not be the same
as the one from Eq. @eq:convection-diffusion-weak. To get rid of this
inconsistency, the anisotropic diffusion can be introduced via a
residual term

$ integral_Omega eps grad u dot.op grad w dif bx + integral_Omega (bold(v) dot.op grad u) w dif bx\
+ sum_(K in cal(M)) delta_K integral_K (- eps Delta u + bold(v) dot.op grad u - f) (bold(v) dot.op grad w) = integral_Omega f (bx) w dif bx $
the added term will be zero for the exact solution (strong PDE) and the
anisotropic diffusion is still here. The control parameter is usually
chosen according to
$ delta_K = cases(delim: "{", eps^(- 1) h_K^2 & quad upright("if ") norm(bold(v))_(K , oo) h_K lt.eq  2 eps, h_K & quad upright("if ") norm(bold(v))_(K , oo) h_K >  2 eps) $
With this, the $Order(h_(cal(M))^2)$ convergence of
$norm(u-u_h)_(L^2 (Omega))$ for $h$-refinement is preserved, while upwind
quadrature only achieves $Order(h_(cal(M)))$ convergence.

== Discretization of Time-Dependent (Transient) Convection-Diffusion IBVPs
<sub:discrete-time-dependent-convection-diffusion>
Now we will take a look at how time-dependent convection-diffusion can
be modeled. Assuming the incompressibility condition and
non-dimensionalizing, @eq:time-dependent-heat becomes
#neq($ frac(partial, partial t) med u - eps Delta u + bold(v \( x ,) t \) dot.op grad u = f quad upright("in ") Omega $) <eq:transient_conv_diff>
If we solve this with the method of lines without upwind quadrature, oscillations occur. But with upwinding, damping is observed, which is  also wrong. 
Therefore, we need a different method. Of course the limit of
$eps arrow.r 0$ again poses a problem. So let's first look at the
pure transport problem, which we get by setting $eps = 0$:
$ frac(partial, partial t) med u + bold(v \( x ,) t \) dot.op grad u = f quad upright("in ") Omega $
Its solution is given by the #emph[Method of Characteristics]
#neq($ u (bx , t) = cases(delim: "{", u_0 (bold(x_0)) + integral_0^t f (bold(y) (s) , s) dif s & upright("if ") bold(y) (s) in Omega quad &forall 0 < s < t, g (bold(y) (s_0) , s_0) + integral_(s_0)^t f (bold(y) (s) , s) dif s quad& upright("if ") bold(y) (s_0) in partial Omega and bold(y) (s) in Omega quad &forall s_0 < s < t) $) <eq:moc>
where $bold(y) (t)$ is defined as the solution of
$dot(bold(y)) (t) = bold(v \( y) (t) , t \), $ $u_0$ the initial
condition and $g$ the Dirichlet boundary conditions on the inflow
boundary. Unfortunately, this only works for the pure transport problem.
For $eps > 0$ we need an other method.

#strong[Splitting Methods]

Given a general ODE whose right-hand side is
the sum of two functions
$ dot(bold(y)) = bold(g) (t , bold(y)) + bold(r) (t , bold(y)) $ The
Strang Splitting single step method provides a method to solve this.

#mybox("Strang Splitting")[
  Compute $bold(y)^((j + 1))$ given
  $bold(y)^((j))$ according to
  $  tilde(bold(y)) &= bold(z) (t_j + 1 / 2 tau) , #h(2em) &&upright(" where ") bold(z) (t) &&upright(" solves ") dot(bold(z)) &&= bold(g) (t , bold(z)) , bold(z) (t_j) = bold(y)^((j))\
  hat(bold(y)) &= bold(w) (t_(j + 1)) , #h(2em) &&upright(" where ") bold(w) (t) &&upright(" solves ") dot(bold(z)) &&= bold(r) (t , bold(w)) , bold(w) (t_j) = tilde(bold(y))\
  bold(y)^((j + 1)) &= bold(z) (t_(j + 1)) , #h(2em) &&upright(" where ") bold(z) (t) &&upright(" solves ") dot(bold(z)) &&= bold(g) (t , bold(z)) , bold(z) (t_j + 1 / 2 tau) = hat(bold(y)) $
  and $t_(j + 1) = t_j + tau$
]
#v(-1cm)
#theorem(number: "10.3.3.5", "Order of Strang splitting single step method")[
  If the IVP in each sub-step is solved exactly or with a 2nd-order time stepping method, the Strang splitting method is of second order.
]
We can now apply this to Eq. @eq:transient_conv_diff.
$ frac(partial, partial t) &u &= quad& eps  Delta u & quad f - bold(v) & dot.op grad u\
  &arrow.t.b & quad & med med arrow.t.b & quad & med med arrow.t.b\
  & dot(bold(y)) &= quad bold(&g (y)) & &bold(r(y)) $ This amounts to once
solving pure diffusion
$ frac(partial t, partial z) - eps Delta z = 0 $ and once pure
transport
$ frac(partial t, partial w) + bold(v) dot.op grad u = f $ To
solve the pure transport problem, we have seen the method of
characteristics @eq:moc. However, it requires integration along
streamlines. We can approximate these integrals by "following" particles along their path defined by the velocity field. This *particle method* works as follows:

+ Pick suitable interpolation nodes $bold(p)_i$, the initial particle
  positions

+ Solve initial value problems
  $ bold(dot(y)) (t) = bold(y)(bold(v)(t), t) quad , quad bold(y) (0) = bold(p)_i $
  with suitable single step methods

+ Reconstruct the approximation. With the composite midpoint rule
  $ u_h^((j)) (bold(p)_i^((j))) = u_0 (bold(p)_i) + tau sum_(l = 1)^(j - 1) f (1 / 2 (bold(p)_i^l + bold(p)_i^(l - 1)) , 1 / 2 (t_l + t_(l - 1))) $

Basically, we have an interpolation problem with nodes that change over time.
At each step, we need to add particles at the inflow boundary and remove
ones that leave the domain. This means we need to re-mesh in each step,
 i.e., create a new triangular mesh with the advected nodes/particles.

#counter(heading).update((10,3,3))
=== Semi-Lagrangian Method
<sub:semi-lagrangian>

#tip-box("")[
  Check out this #weblink("https://youtu.be/kvBRFxRIJuY", "video") for an intuitive explanation of the method in 1D.
]

A method which relies on a _fixed_ mesh is the Semi-Lagrangian method.
#definition(number: "10.3.4.2", "Material derivative", ..unimportant)[
  Given a velocity field $bold(v)$, the material derivative of a function
  $f$ is given by
  $ frac(D f, D bold(v)) (bx , t_0) = lim_(tau arrow.r 0) frac(f (bx , t_0) - f (Phi^(t_0 , t_0 - tau) bx , t_0 - tau), tau) $
]

By the chain rule we find
$ frac(D f, D bold(v)) (bx , t) = grad f (bx , t) dot.op bold(v \( x) , t \) + frac(partial t, partial f) (bx , t) $

Hence the transient convection-diffusion equation
@eq:transient_conv_diff can be rewritten as
$ frac(D u, D bold(v)) - eps Delta u = f quad upright(" in ") Omega $
By using a backwards difference of the material derivative, we get a
semi-discretization
$ frac(u^((j)) (bx) - u^((j - 1)) (Phi^(t_j , t_j - tau) bx), tau) - eps Delta u^((j)) = f (bx , t_j) quad upright(" in ") Omega $
with additional initial conditions for $t = t_j$. On this
semi-discretization the standard Galerkin method can be applied.
$ integral_Omega frac(u^((j)) (bx) - u^((j - 1)) (Phi^(t_j , t_j - tau) bx), tau) w(bx) dif bx + eps integral_Omega grad u^((j)) dot.op grad w dif bx = integral_Omega f (bx , t_j) w(bx) dif bx $
Unfortunately this cannot be implemented as is, because the function
$u^((j - 1)) (Phi^(t_j , t_j - tau) bx)$ is not smooth in $cal(M)$
and is hence not a finite element function on $cal(M)$. To get around
this, simply replace it by its linear interpolant
$I_1 (u^((j - 1)) circle.stroked.tiny Phi^(t_j , t_j - tau))$ and
replace $Phi^(t_j , t_j - tau) bx$ by
$bx - tau bold(v)(bx, t_j \)$ (explicit Euler).
$ integral_Omega frac(u^((j)) (bx) - I_1 lr({u^((j - 1)) (bold(s) - tau bold(v)(bold(s) , t_j))}, size: #110%) (bx), tau) w(bx) dif bx + eps integral_Omega grad u^((j)) dot.op grad w dif bx = integral_Omega f (bx , t_j) w(bx) dif bx $
Here, $bold(s)$ is the variable on which the interpolation operator $I_1$ acts (in the lecture notes, this variable is written as a dot "$med dot.op med$"). The above equation still looks quite complicated, but can be implemented.
