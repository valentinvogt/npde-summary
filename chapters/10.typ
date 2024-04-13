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
$ bold(j (x)) = - kappa grad med u (bx) + bold(v (x)) rho u (bx) $
with $kappa$ the heat conductivity and $rho$ the volumetric heat
capacity. We already know the first part, called diffusive heat flux, from the
heat equation. The second part is called convective heat flux. With this
new flux, the standard PDE becomes
#neq($ - div (kappa grad med u) + div (rho bold(v (x)) u) = f quad upright("in ") Omega $) <eq:convection-diffusion-strong>
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
$ - div (kappa grad med u) + bold(v (x)) dot.op grad rho u = f quad upright("in ") Omega $
We can also look at the time-dependent heat flow, which is similar to
the heat equation @eq:wave-variational
$ frac(partial, partial t) med (rho u) - div (kappa grad med u) + div (rho med bold(v) (bx ,t) med bu) = f quad upright("in ") Omega $
We will see later on how this can be solved.

#pagebreak()
== Stationary Convection-Diffusion Problems
<sub:stationary-convection-diffusion>
Here we will focus on the convection diffusion equation
@eq:convection-diffusion-strong with constant $kappa$, $rho$, incompressible flow
$bold(v)$ and zero Dirichlet boundary conditions. Hence
$ - kappa Delta u + rho bold(v (x)) dot.op grad u = f quad upright("in ") Omega , #h(2em) u = 0 quad upright("on ") partial Omega $
Non-dimensionalizing the problem results in
$ - epsilon.alt Delta u + bold(v (x)) dot.op grad u = f quad upright("in ") Omega , #h(2em) u = 0 quad upright("on ") partial Omega $
with $norm(bold(v))_(L^oo (Omega)) = 1$. This results in the following
variational form
#neq($ epsilon.alt integral_Omega grad u dot.op grad v dif bx + integral_Omega (bold(v) dot.op grad u) v dif bx = integral_Omega f (bx) v dif x $) <eq:convection-diffusion-weak>
with the left-hand side the bilinear form $a (u , v)$. However, $a$ is
not symmetric. This also means that it does not induce an energy norm.
However, it is still positive definite (see lecture document).

#strong[Singular perturbation] A boundary value problem depending on a
parameter $epsilon.alt$ is called singularly perturbed, if the limit
problem for $epsilon.alt arrow.r epsilon.alt_0$ is not compatible with
the boundary conditions.

For $epsilon.alt = 0$ the above PDE is singular perturbed. It cannot
satisfy Dirichlet boundary conditions on the outflow part of the
boundary.
$Gamma_(upright("out")) = bx in partial Omega : bold(v (x)) dot.op bold(n (x)) > 0$,
similarly
$Gamma_(upright("in")) = bx in partial Omega : bold(v (x)) dot.op bold(n (x)) < 0$

#strong[Upwinding]

When trying to solve Eq. @eq:convection-diffusion-weak with
the Galerkin approach, when $epsilon.alt$ is very close to 0, one can
observe huge oscillations in the solution, which is not correct. It
comes from the fact that the Galerkin matrix becomes close to singular.
So our goal is to get a robust method that can solve Eq.
@eq:convection-diffusion-weak no matter the $epsilon.alt$.

Consider again Eq. @eq:convection-diffusion-weak but in $d = 1$ and with zero
boundary conditions
$ epsilon.alt integral_0^1 frac(partial u, partial x) frac(partial v, partial x) dif x + integral_0^1 frac(partial u, partial x) v dif x = integral_0^1 f (x) v dif x $
To calculate the Galerkin matrix for an equidistant mesh with $M$ cells,
we use the global composite trapezoidal rule for the convective term
$ integral_0^1 psi (x) dif x = h sum_(j = 0)^M psi (j h) $
Hence the convective term of the bilinear form will be approximated by
$ integral_0^1 frac(partial u_h, partial x) v_h dif x approx h sum_(j = 1)^(M - 1) frac(partial u_h, partial x) med (j h) v_h (j h) $
But $frac(partial u_h, partial x) med (j h)$ is not valid, as it's
discontinuous at the nodes for $u_h in cal(S)_(1 , 0)^0$. However,
convection transports the information in the direction of $bold(v)$ ($= 1$
in our case). Hence use
$ frac(partial u_h, partial x) med (j h) = lim_(delta arrow.r 0) frac(partial u_h, partial x) med (j h - delta bold(v)) = lr(frac(partial u_h, partial x) mid(|)) _( \] x_(j - 1) , x_j \[) $ 
And generalized in more dimensions
$ bold(v (p)) dot.op grad u_h (bold(p)) = lim_(delta arrow.r 0) bold(v (p)) dot.op grad u_h (bold(p) - delta bold(v (p))) $

#strong[Streamline diffusion]

A totally different idea to fix the
problem of $epsilon.alt arrow.r 0$ is to add some $h$-dependent
diffusion. I.e., replace $epsilon.alt arrow.l epsilon.alt + c (h)$ with
$c (h) > 0$. However, there is smearing in the internal layers. But as
the solution is smooth along the direction of $bold(v)$, so adding
diffusion along the velocity should not do any harm.

The method of #emph[Anisotropic diffusion] is born. On cell $K$ replace
$epsilon.alt arrow.l epsilon.alt bold(I) + delta_K bold(v)_K bold(v)_K^top$
with $bold(v)_K$ the local velocity, i.e. obtained by averaging and
$delta_K > 0$ some controlling parameter. Resulting in
$ integral_Omega (epsilon.alt bold(I) + delta_K bold(v)_K bold(v)_K^top) grad u dot.op grad v dif bx + integral_Omega (bold(v) dot.op grad u) v dif bx = integral_Omega f (bx) v dif bx $
However this affects the solution $u$, such that it will not be the same
as the one from Eq. @eq:convection-diffusion-weak. To get rid of this
inconsistency, the anisotropic diffusion can be introduced via a
residual term
$ integral_Omega epsilon.alt grad u dot.op grad v dif bx + integral_Omega (bold(v) dot.op grad u) v dif bx\
+ sum_(K in cal(M)) delta_K integral_K (- epsilon.alt Delta + bold(v) dot.op grad u - f) bold(v) dot.op grad v = integral_Omega f (bx) v dif bx $
the added term will be zero for the exact solution (strong PDE) and the
anisotropic diffusion is still here. The control parameter is usually
chosen according to
$ delta_K = cases(delim: "{", epsilon.alt^(- 1) h_K^2 & upright("if ") norm(bold(v))_(K , oo) h_K lt.eq  2 epsilon.alt, h & upright("if ") norm(bold(v))_(K , oo) h_K >  2 epsilon.alt) $
With this, the $cal(O) (h_(cal(M))^2)$ convergence of
$norm(u-u_h)_(L^2 (Omega))$ for $h$ -efinement is preserved, while upwind
quadrature only achieves $cal(O) (h_(cal(M)))$ convergence.

== Discretization of Time-Dependent (Transient) Convection-Diffusion IBVPs
<sub:discrete-time-dependent-convection-diffusion>
Now we will take a look at how time dependent convection diffusion can
be modeled. Assuming the incompressibility condition and
non-dimensionalizing, Eq. @eq:heat-integral-form becomes
#neq($ frac(partial, partial t) med u - epsilon.alt Delta u + bold(v \( x ,) t \) dot.op grad u = f quad upright("in ") Omega $) <eq:transient_conv_diff>
Up on inspecting the solution obtained with method of lines, one
observes that without upwind quadrature, oscillations occur. However,
with upwind damping is observed, which is wrong. Hence other methods of
solving have to be explored. Of course the limit of
$epsilon.alt arrow.r 0$ again poses a problem. So lets first look at the
pure transport problem
$ frac(partial, partial t) med u + bold(v \( x ,) t \) dot.op grad u = f quad upright("in ") Omega $
Its solution is given by the #emph[Method of Characteristics]
#neq($ u (bx , t) = cases(delim: "{", u_0 (bold(x_0)) + integral_0^t f (bold(y) (s) , s) dif s & upright("if ") bold(y) (s) in Omega quad &forall 0 < s < t, g (bold(y) (s_0) , s_0) + integral_(s_0)^t f (bold(y) (s) , s) dif s quad& upright("if ") bold(y) (s_0) in partial Omega and bold(y) (s) in Omega quad &forall s_0 < s < t) $) <eq:moc>
where $ dot(bold(y)) (t) = bold(v \( y) (t) , t \) $ $u_0$ the initial
condition and $g$ the Dirichlet boundary conditions on the inflow
boundary. Unfortunately, this only works for the pure transport problem.
For $epsilon.alt > 0$ we need an other method,

#strong[Splitting Methods]

Given a general ODE whose right-hand side is
the sum of two functions
$ dot(bold(y)) = bold(g) (t , bold(y)) + bold(r) (t , bold(y)) $ The
Strang Splitting single step method provides a method to solve this.

#mybox("Strang Splitting")[
  Compute $bold(y)^((j + 1))$ given
  $bold(y)^((j))$ according to
  $  tilde(bold(y)) = bold(z) (t_j + 1 / 2 tau) , #h(2em) upright(" where ") bold(z) (t) upright(" solves ") dot(bold(z)) = bold(g) (t , bold(z)) , bold(z) (t_j) = bold(y)^((j))\
  hat(bold(y)) = bold(w) (t_(j + 1)) , #h(2em) upright(" where ") bold(w) (t) upright(" solves ") dot(bold(z)) = bold(r) (t , bold(w)) , bold(w) (t_j) = tilde(bold(y))\
  bold(y)^((j + 1)) = bold(z) (t_(j + 1)) , #h(2em) upright(" where ") bold(z) (t) upright(" solves ") dot(bold(z)) = bold(g) (t , bold(z)) , bold(z) (t_j + 1 / 2) = hat(bold(y)) $
  and $t_(j + 1) = t_j + tau$
]
#v(-1cm)
#theorem(number: "10.3.3.5", "Order of Strang splitting single step method")[
  Assuming exact or second order accuracy solution of the initial value
  problems of the sub-steps, the Strang splitting method is of
  second order.
]
We can now apply this to Eq. @eq:transient_conv_diff.
$ A && B quad && C E \
 M quad && N && P $
$ frac(partial, partial t) &u &= quad& epsilon.alt  Delta u & quad f - bold(v) & dot.op grad u\
  &arrow.t.b & quad & med med arrow.t.b & quad & med med arrow.t.b\
  & dot(bold(y)) &= quad bold(&g (y)) & &bold(r(y)) $ This amount to once
// $ frac(partial, partial t) &u quad &= & quad epsilon.alt Delta u & quad f - bold(v) dot.op grad u\
// arrow.t.b &  & arrow.t.b & arrow.t.b\
// dot(bold(&y)) & =  bold(g (y)) & bold(r (y)) $ This amount to once
solving pure diffusion
$ frac(partial t, partial z) - epsilon.alt Delta z = 0 $ and once pure
transport
$ frac(partial t, partial w) + bold(v) dot.op grad u = f $ To
solve the pure transport problem, we have seen the method of
characteristics @eq:moc. However, it requires integration along
streamlines. One idea is to solve it with the particle method.

+ Pick suitable interpolation nodes $bold(p)_i$, the initial particle
  positions

+ Solve initial value problems
  $ bold(dot(y)) (t) = bold(v \( y ,) t \) quad , quad bold(y) (0) = bold(p)_i $
  with suitable single step methods

+ Reconstruct the approximation. With the composite midpoint rule
  $ u_h^((j)) (bold(p)_i^((j))) = u_0 (bold(p)_i) + tau sum_(l = 1)^(j - 1) f (1 / 2 (bold(p)_i^l + bold(p)_i^(l - 1)) , 1 / 2 (t_l + t_(l - 1))) $

But the interpolation nodes change over time and care needs to be taken,
to add particles each step at the inflow boundary and remove ones, which
leave the domain. Because of the movement of the nodes and potential
creation and deletion, each step we need to re-mesh, i.e., create a new
triangular mesh with the advected nodes/particles.

#strong[Semi Lagrangian]

Another method, which relies on a fixed mesh, is the semi Lagrangian method.
#definition(number: "10.3.4.2", "Material derivative")[
  Given a velocity field $bold(v)$, the material derivative of a function
  $f$ is given by
  $ frac(D f, D bold(v)) (bx , t_0) = lim_(tau arrow.r 0) frac(f (bx , t_0) - f (Phi^(t_0 , - tau) bx , t_0 - tau), tau) $
]

By the chain rule we find
$ frac(D f, D bold(v)) (bx , t) = grad f (bx , t) dot.op bold(v \( x) , t \) + frac(partial t, partial f) (bx , t) $
Hence the transient convection diffusion Eq.
@eq:transient_conv_diff can be rewritten as
$ frac(D u, D bold(v)) - epsilon.alt Delta u = f quad upright(" in ") Omega $
By using a backwards difference of the material derivative, we get a
semi-discretization
$ frac(u^((j)) (bx) - u^((j - 1)) (Phi^(t_j , t_j - tau) bx), tau) - epsilon.alt Delta u^((j)) = f (bx , t_j) quad upright(" in ") Omega $
with additional initial conditions for $t = t_j$. On this
semi-discretization the standard Galerkin method can be applied.
$ integral_Omega frac(u^((j)) (bx) - u^((j - 1)) (Phi^(t_j , t_j - tau) bx), tau) v dif bx + epsilon.alt integral_Omega grad u^((j)) dot.op grad v dif bx = integral_Omega f (bx , t_j) v dif bx $
Unfortunately this cannot be implemented as is, because the function
$u^((j - 1)) (Phi^(t_j , t_j - tau) bx)$ is not smooth in $cal(M)$
and is hence not a finite element function on $cal(M)$. To get around
this, simply replace it by its linear interpolant
$I_1 (u^((j - 1)) circle.stroked.tiny Phi^(t_j , t_j - tau))$ and
replace $Phi^(t_j , t_j - tau) bx$ by
$bx - tau bold(v \( x ,) t_j \)$ (explicit Euler).
$ integral_Omega frac(u^((j)) (bx) - I_i (u^((j - 1)) (dot.op - tau bold(v \( dot.op ,) t_j)) \) (bx), tau) v dif bx + epsilon.alt integral_Omega grad u^((j)) dot.op grad v dif bx = integral_Omega f (bx , t_j) v dif bx $
This can be implemented now.
