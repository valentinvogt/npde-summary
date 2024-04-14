#import "../src/setup.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#show: thmrules

= Finite Element Method
<ch:finite-element-method>
#counter(heading).step(level: 2)
== Galerkin Discretization
<sub:galerkin-discretization>
The idea is to replace the infinite function space $V_0$ by a finite subspace $V_(0 , h) subset V_0$

#theorem(number: "2.2.1.5", title: "Theorem", "Unique solution of discrete variational problems")[
  If the bilinear form $a : V_0 times V_0 arrow.r bb(R)$ is symmetric and positive definite and the linear form $ell : V_0 arrow.r bb(R)$ is continuous, then the discrete variational problem:
  #neq($ u_h in V_(0 , h) : a (u_h , v_h) = ell (v_h) , quad forall v_h in V_(0 , h) $)<eq:disc-var-prob>
  has a unique #emph[Galerkin] solution $u_h in V_(0 , h)$ satisfying the energy estimate
  $ norm(u)_a lt.eq sup_(v_h in V_(0 , h)) frac(lr(|ell (v_h)|), norm(v_h)_a) $
] <thm:galerkin-existence-and-uniqueness>

Recall the definition of a basis (here, superscripts are indices and not to be confused with exponents): ${ b^1 , dots.h , b^N } subset V$ is
a basis, if for every $v in V$ there are unique coefficients $mu_l$ such
that $v = sum_(i = 1)^N mu_i b^i$ and $N$ agrees with the dimension of
$V$. Now we can expand $u_h = mu_1 b^1 + dots.h.c + mu_N b^N$ and our
goal is to find the coefficients $mu_i$.

#mybox("Galerkin Discretization")[
  Linear discrete variational problem @eq:disc-var-prob $xarrow(sym: -->, "Choice of basis" frak(B)_h)$
  Linear system of equations $bold(A arrow(mu) = arrow(phi))$
  $ bold("Galerkin Matrix") : & quad bA = [a (b^k , b^j)]_(j , k = 1)^N in bb(R)^(N , N)\
  bold("RHS vector") : & quad bold(arrow(phi)) = [ell (b^j)]_(j = 1)^N in bb(R)^N\
  bold("coefficient vector") : & quad bold(arrow(mu)) = [mu_1 , dots.h , mu_N]^top in bb(R)^N $
]

Note that $bA_(j , k) = a (b_h^k , b_h^j) eq.not a (b_h^j , b_h^k)$
if $a$ is not symmetric (note the order $k,j$). Of course the bilinear form $a$ determines
some properties of the Galerkin matrix. If $a$ is symmetric and/or
positive definite, the Galerkin matrix $bA$ will have the same
properties.

The choice of $V_(0 , h)$ alone determines the quality of the solution
$u_h$. While mathematically the choice of basis $frak(B)_h$ does not
matter, for solving the equation numerically, the choice is crucial as
the basis determines how stable and efficiently the solution can be
computed, as it determines for example the sparsity of $bA$.

#tip([Computing the energy norm in code])[
  Sometimes, problems ask you to compute $norm(u_h)_a$, i.e., the energy norm of a discrete solution. 
  // $|u_h|_a = a(u_h,u_h) = a\left(\sum_{i=0}^N \mu_i \cdot b_i^h, \sum_{i=0}^N \mu_i \cdot b_i^h\right) $
  $ norm(u_h)_a = sqrt(a(u_h,u_h)) = a(sum_(i=0)^N mu_i b_i^h, sum_(j=0)^N mu_j b_j^h) = sum_(i=0)^N sum_(j=0)^N mu_i mu_j a(b_i^h,b_j^h) = bold(arrow(mu))^top bA bold(arrow(mu)) $
  So it can be computed as `sqrt(mu.transpose() * A * mu)`.
]

#pagebreak()
== Linear FEM in 1D
<sub:linfem1d>
In FEM, the goal is to approximate $u$ by piecewise polynomial functions.

#mybox("Mesh in one dimension")[
  Let $Omega = [a , b]$, we equip it with $M + 1$ #strong[nodes] resulting in the set of nodes:
  $ cal(V (M)) = { a = x_0 < x 1 < dots.h.c < x_M = b } $
  The nodes define intervals, which build up the mesh:
  $ cal(M) = { \] x_(j - 1) , x_j \[ : 1 lt.eq j lt.eq M } $
  The intervals $[x_(j - 1) , x_j]$ are the #strong[cells] of the mesh.
  We define local cell size $h_j = lr(|x_j - x_(j - 1)|)$ and global mesh width $h_(cal(M)) = max_j h_j$.
]

A simple space for continuous, $cal(M)$-piecewise polynomial functions
in $H_0^1 (\] a , b \[)$ consists of linear functions on each cell:
#neq($ V_(0 , h) = S_(1 , 0)^0 (cal(M)) = {v in C^0 ([a , b]) : eval(v)_[x_(i - 1),x_i ] upright("is linear") , i = 1 , dots.h , M , v (a) = v (b) = 0} $)<eq:linfem1d-space>
$ arrow.r N = upright("dim") S_(1 , 0)^0 (cal(M)) = M - 1 $ The
0-superscript stands for global $C^0$ of the functions. The 1-subscript
denotes local degree 1 polynomial and the 0-subscript denotes value 0 on the
boundary. The $cal(S)$ stands for $cal(S)$calar functions.

Common basis functions are the 1D tent functions:
#neq($ b_h^j (x) = cases(delim: "{", (x - x_(j - 1)) \/ h_j & upright("if ") x_(j - 1) lt.eq x lt.eq x_j, (x_j - x) \/ h_(j + 1) & upright("if ") x_j lt.eq x lt.eq x_(j + 1), 0 & "else") $)
$ arrow.r b_h^j (x_i) = delta_(i j) $ A basis satisfying condition (25)
is called a #strong[cardinal] basis. Another key property of tent
functions is that their support just comprises two adjacent cells:
$ "supp"(b_h^j) = [x_(j - 1) , x_(j + 1)] $ The support of a function $f : Omega arrow.r bb(R)$ is defined as $overline({ x in Omega : f (x) eq.not 0 })$, the set of inputs for which the function is nonzero.
Since polynomials are easy to differentiate and integrate, computing (bi)linear forms $a$ and $ell$ for them is quite easy.

#pagebreak()
== Linear FEM in 2D
<sub:linfem2d>
#mybox("Mesh in two dimension")[
  Meshes in 2D rely on #strong[triangulations];. A #strong[triangulation] $cal(M)$ of $Omega$ satisfies:

  - $cal(M) = { K_i }$, where $K_i$ are open triangles

  - $i eq.not j arrow.r K_i sect K_j = diameter$

  - $union.big_(i = 1)^M overline(K_i) = overline(Omega)$

  - $i eq.not j arrow.r overline(K_i) sect overline(K_j)$ is either
    $diameter$, an edge from both triangles or a vertex from both

  Again, the vertices are called #strong[nodes] and the triangles are the
  #strong[cells];.
]


#grid(
  columns: (0.75fr, 0.2fr),
  "This definition does not allow for hanging nodes because of point 4. Hanging nodes are those which lie on the edge of a triangle:", 
  image("../images/hanging_node.png", width: 80%),
)
#v(-0.9cm)
We define a space of piecewise-linear functions analogously to @eq:linfem1d-space:
#neq($ V_(0 , h) = S_1^0 (cal(M)) = {v in C^0 (overline(Omega)) : v_K (bx) = alpha_K + bold(beta)_K dot.op bx , alpha_K in bb(R) , bold(beta)_K in bb(R)^2 , bx in K} $)<eq:linfem2d-space>
$ arrow.r upright("dim") S_1^0 (cal(M)) = \# cal(V (M)) $ And
$S_(1 , 0)^0 (cal(M))$ additionally requires functions to be zero
on $partial Omega$, with
$ upright("dim") S_(1 , 0)^0 (cal(M)) = \# { x in cal(V (M)) : in.not partial Omega } } $
Similarly, the 1D tent functions can be extended to 2D by requiring the
cardinal property. This property is already enough since three points
fully define a plane — in other words, knowing the values in three
points (the vertices of a triangle) fully defines a linear function.
Cardinal bases will produce sparse Galerkin matrices, as the support of a basis function only covers the triangles adjacent to its node, which can only interact with neighboring basis functions.

#mybox("Computation of Galerkin Matrix")[
  Often bilinear forms are defined by integration over the whole domain. But we have seen that the support of basis functions is only local. We can exploit this by performing integration only over the cells.
  #neq($ bA_(i j) = a (b_h^j , b_h^i) = sum_(K in "supp"(b_h^j) sect "supp"(b_h^i)) eval(a)_K (b_h^j , b_h^i) $)<eq:galerkin_matrix_entry>
  where $eval(a)_K$ is the local bilinear form over cell $K$. 
]
#v(-1.1cm)
#mybox("Cell oriented assembly")[
  To further take advantage of @eq:galerkin_matrix_entry, cell oriented assembly can be performed. Go through all cells and compute $eval(a)_K  (b_h^j , b_h^i)$ for all pairs of basis functions associated with cell $K$ (element matrix) and add the values to the entry at $(i , j)$ of $bA$.
]
The same procedure can be applied to calculating the right hand side
vector $phi$, just that only one basis function is involved as the RHS comes from a linear functional.

== Building Blocks of General Finite Element Methods
<sub:fem-building-blocks>
The first building block are meshes, see @sub:linfem1d and @sub:linfem2d.
Next we need to choose a space of functions.

#definition(number: "2.5.2.2", "Multivariate Polynomials")[
  Space of d-variate (taking inputs in $RR^d$) polynomials with total degree $p$:
  $ cal(P)_p (bb(R)^d) = {bx in bb(R)^d arrow.r sum_(bold(alpha) in bb(N)_0^d , lr(|alpha|) lt.eq p) c_alpha bx^alpha , c_alpha in bb(R)} $
  with
  $alpha = (alpha_1 , dots.h.c , alpha_d) , bx^alpha = x_1^(alpha_1) dot.op dots.h.c dot.op x_d^(alpha_d)$
  and $lr(|alpha|) = alpha_1 + dots.h.c + alpha_d$
] <concept:multi-index>

As an example,
$cal(P)_2 (bb(R)^2) = "Span"{ 1 , x_1 , x_2 , x_1^2 , x_2^2 , x_1 x_2 }$

#lemma(number: "2.5.2.5", "Dimension of spaces of Polynomials")[
  $ upright("dim") cal(P)_p (bb(R)^d) = binom(d + p, p) , quad p in bb(N)_0 , d in bb(N) $
  which in the limit of $p arrow.r oo$ behaves like $Order(p^d)$.
]

#definition(number: "2.5.2.7", "Tensor Product Polynomials")[
  The space of tensor product polynomials of degree $p$ in each coordinate is
  $ cal(Q)_p (bb(R)^d) = {bx in bb(R)^d arrow.r sum_(l_1 = 0)^p dots.h.c sum_(l_d = 0)^p c_(l_1 , dots.h , l_d) x_1^(l_1) dot.op dots.h.c dot.op x_d^(l_d) , quad "where" c_(l_1 , dots.h , l_d) in bb(R)}\
  = "Span"{bx arrow.r p_1 (x_1) dot.op dots.h.c dot.op p_d (x_d) , p_i in cal(P)_p (bb(R))} $

]
As an example,
$cal(Q)_2 (bb(R)^2) = "Span" { 1 , x_1 , x_2 , x_1 x_2 , x_1^2 , x_1^2 x_2 , x_1^2 x_2^2 , x_1 x_2^2 , x_2^2 }$

#lemma(number: "2.5.2.8", "Dimension of spaces of tensor product polynomials")[
  $ upright("dim") cal(Q)_p (bb(R)^d) = (p + 1)^d , quad p in bb(N)_0 , d in bb(N) $
]
Finally we need #strong[locally] supported basis functions. This basis
$frak(B)_h = b_h^1 , dots.h.c , b_h^N$ should satisfy the following
constraints:

+ $frak(B)_h$ is a basis of $V_h$, hence
  $upright("dim") frak(B)_h = upright("dim") V_h$.

+ Each $b_h^i$ is associated with a single mesh geometry entity
  (cell/edge/face/vertex).

+ Each $b_h^i$ is only locally supported, i.e., only nonzero in adjacent
  cells.

== Lagrangian Finite Element Methods
<sub:lagrangian-fem>
Remember Eqs. @eq:linfem1d-space and @eq:linfem2d-space as two examples of finite element spaces. They are examples of the general Lagrangian FE spaces. First we introduce them for #emph[simplicial] meshes, i.e., those consisting of triangles (2D) or tetrahedra (3D).

#definition(number: "2.6.1.1", "Simplicial Lagrangian finite element spaces")[
  #neq($ cal(S)_p^0 (cal(M)) = {v in C^0 (overline(Omega)) : v_K in cal(P)_p (K) , forall K in cal(M)} $)<eq:simp_lfes>
]
This space is well suited for triangular meshes, as the local dimension
$binom(d + p, p) = binom(2 + p, p)$ (in 2D) is the same as the amount of
vertices in a triangle and its interpolation nodes. The local basis
functions of $cal(S)_1^0$ are the barycentric coordinate functions. In
$cal(S)_2^0$, the local basis functions are linear combinations of
barycentric coordinate functions:
#grid(
  columns: (0.5fr, 0.5fr),
  [
    #h(1cm)
  $ b_K^1 & = (2 lambda_1 - 1) lambda_1 , quad b_K^4 = 4 lambda_1 lambda_2 ,\
b_K^2 & = (2 lambda_2 - 1) lambda_2 , quad b_K^5 = 4 lambda_2 lambda_3 ,\
b_K^3 & = (2 lambda_3 - 1) lambda_3 , quad b_K^6 = 4 lambda_1 lambda_3 $
  ],
  image("../images/triangle_nodes.png", width: 80%),
)
where the local basis functions 1-3 are associated with vertices and 4-6 with edges.

Analogously, the following space is well suited for quadrilaterals, as its local dimension $(p + 1)^2$ is the same as the amount of vertices and interpolation points on quads.

#definition(number: "2.6.2.5", "Tensor product Lagrangian finite element spaces")[
  #neq($ cal(S)_p^0 (cal(M)) = {v in C^0 (overline(Omega)) : eval(v)_K in cal(Q)_p (K) , forall K in cal(M)} $) <eq:tens_lfes>
]

Note that the choice of local polynomial space is the only difference,
$cal(Q)_p (K)$ instead of $cal(P)_p (K)$. Of course these spaces can be mixed: on mixed (so-called hybrid) meshes, we use
@eq:simp_lfes on triangles and @eq:tens_lfes on quadrilaterals.

== Implementation of Finite Element Methods
<sub:implementation-of-fem>
Remember the principle of cell-oriented assembly. The goal is to rely
mostly on local computations. To perform cell-oriented assembly, a map
from local to global indices is needed. In LehrFEM++ this is the job of
the dofhandler (
  #weblink("https://craffael.github.io/lehrfempp/classlf_1_1assemble_1_1_dof_handler.html")[`lf::assemble::DofHandler`];).
It provides the following main methods:

- `NumDofs()`, returns the total number of global basis functions, the
  dimension of the FE space.

- `NumLocalDofs(const lf::mesh::Entity &)`, returns the number of global
  basis functions covering any geometric entity.

- `GlobalDofIndices(const lf::mesh::Entity &)`, returns an array of
  indices of the global basis function covering the given entity. Their order is determined by the local numbering of the basis functions. E.g., the global index of local basis function 0 is given by `GlobalDofIndices(cell)[0]`.

- `NumInteriorDofs(const lf::mesh::Entity &)`, returns the number of
  global basis functions *associated with* the given entity (so for a triangle, the number of basis functions in the interior of the triangle, not those of the edges or vertices).

- `InteriorGlobalDofIndices(const lf::mesh::Entity &)`, similar to
  `GlobalDofIndices` but returns only the indices of the global basis
  functions associated with the given entity.

- `Entity(gdof_idx_t dofn)`, returning the entity associated with the
  global index `dofnum`.

Instead of dimension, in LehrFEM++ the concept of #strong[co-dimension]
is used. Instead of going from a point with dimension 0 to a triangle
with dimension 2, the co-dimension is the other way around. The
highest-dimensional entity has co-dimension 0. This ensures that cells
are always of co-dimension 0.

To assemble the Galerkin matrix,
#weblink("https://craffael.github.io/lehrfempp/group__assemble__matrix__locally.html#ga39b4197203dd4e896bd7073fc033aca3")[`lf::assemble::AssembleMatrixLocally`];
can be used. To use it, we need element matrix providers
(#weblink("https://craffael.github.io/lehrfempp/group__entity__matrix__provider.html")[docs];).
These are constructs which provide the element matrix for given bilinear
forms. Some common bilinear forms are already implemented.

- $integral_K alpha (bx) med grad med u dot.op grad med v dif bx$
  is implemented in
  #weblink("https://craffael.github.io/lehrfempp/classlf_1_1fe_1_1_diffusion_element_matrix_provider.html")[`lf::fe::DiffusionElementMatrixProvider`]

- $integral_K gamma (bx) med u med v dif bx$ is
  implemented in
  #weblink("https://craffael.github.io/lehrfempp/classlf_1_1fe_1_1_mass_element_matrix_provider.html")[`lf::fe::MassElementMatrixProvider`]

- $integral_e gamma (bx) med u med v dif S$ is
  implemented in
  #weblink("https://craffael.github.io/lehrfempp/classlf_1_1fe_1_1_mass_edge_matrix_provider.html")[`lf::fe::MassEdgeMatrixProvider`];.
  Note the integration over an edge and not a cell.

- $integral_K alpha (bx) med grad med u dot.op grad med v dif bx med med + integral_K gamma (bx) med u med v dif bx$
  combined is implemented in \
  #weblink("https://craffael.github.io/lehrfempp/classlf_1_1uscalfe_1_1_reaction_diffusion_element_matrix_provider.html")[`lf::uscalfe::ReactionDiffusionElementMatrixProvider`]

- $integral_K f (bx) v dif bx$ is implemented by
  #weblink("https://craffael.github.io/lehrfempp/classlf_1_1fe_1_1_scalar_load_element_vector_provider.html")[`lf::fe::ScalarLoadElementVectorProvider`]

- $integral_e f (bx) v  dif S$ is implemented by
  #weblink("https://craffael.github.io/lehrfempp/classlf_1_1fe_1_1_scalar_load_edge_vector_provider.html")[`lf::fe::ScalarLoadEdgeVectorProvider`];.
  Note again the integration over an edge.

Note that the last two are actually element vector providers.

#lemma(number: "2.7.5.5", "Integration of powers of barycentric coordinate functions")[
  For a $d$-simplex $K$ (line in 1D, triangle in 2D, tetrahedron in 3D) with barycentric coordinate functions
  $lambda_1 , dots.h , lambda_(d + 1)$
  #neq($ integral_K lambda_1^(alpha_1) dot.op dots.h.c dot.op lambda_(d + 1)^(alpha_(d + 1)) dif bx = d ! lr(|K|) frac(alpha_1 ! dot.op dots.h.c dot.op alpha_(d + 1) !, (alpha_1 + dots.h.c + alpha_(d + 1) + d) !) , quad alpha_i in bb(N), $)
  where $|K|$ is the volume of the simplex.
]

#strong[Quadrature rules]
$ integral_K f (x) dif x approx sum_(l = 1)^(P_K) w_l^K f (zeta_l^K) , quad w_l^K arrow.r upright("weights") , zeta_l^K arrow.r upright("(quadrature) nodes") $
#subtle-box()[
  *Order* of a quadrature rule: a quad rule is of order $q$ if

  - for a simplex $K$, it is exact for all polynomials
    $f in cal(P)_(p - 1) (bb(R)^d)$

  - for a tensor product element $K$, it is exact for all polynomials
    $f in cal(Q)_(p - 1) (bb(R)^d)$
]

#strong[Essential Boundary Conditions] Remember from @sub:essential-and-natural-boundary-conditions that essential
boundary conditions are Dirichlet boundary conditions, i.e., $u = g$ on
$partial Omega$, and can be solved with the offset function trick. This
trick can also be used in FEM. \
Assume that the basis functions are sorted such that all interior ones
come first, followed by the ones on the boundary (we are free to choose
the order of basis functions). Then, it is possible to write $A$ as
follows:
$ bA = mat(delim: "[", bA_0, bA_(0 partial); bA_(0 partial)^top, bA_(partial partial)) $
where $bA_0$ is the Galerkin matrix for $cal(S)_(p , 0)^0 (cal(M))$,
containing the interactions among interior basis functions.
$(bA_(0 partial))_(i j) = a (b_h^j , b_h^i)$, where $b_h^j$ belongs to
the boundary and $b_h^i$ to the interior, so $bA_(0 partial)$ contains
the interactions of interior with boundary functions. Similarly,
$bA_(partial partial)$ consists only of entries calculated from basis
functions of the boundary. Then we want to solve:
$ mat(delim: "[", bA_0, bA_(0 partial); bA_(0 partial)^top, bA_(partial partial)) mat(delim: "[", bold(arrow(mu))_0; bold(arrow(mu))_partial) = mat(delim: "[", bold(arrow(phi)); bold(arrow(phi))_partial) $
where $mu_partial$ are the coefficients of the basis expansion of $g$ on
the boundary, which we know since $g$ is given. We only need to solve
for $mu_0$ which results in
$ bA_0 bold(arrow(mu))_0 = bold(arrow(phi)) - bA_(0 partial) bold(arrow(mu))_partial $

This can be done in LehrFEM++ with
#weblink("https://craffael.github.io/lehrfempp/namespacelf_1_1assemble.html#a4fba0f99e10227530fcb990ddda7b305")[lf::assemble::FixFlaggedSolutionComponents]
or \
#weblink("https://craffael.github.io/lehrfempp/namespacelf_1_1assemble.html#ad8de42b7c7e79eeba5704e43a5b4d67f")[lf::assemble::FixFlaggedSolutionCompAlt];.
Both modify the matrix $bA$ and RHS vector $bold(arrow(b))$ such that $bold(A arrow(mu) = arrow(b))$ has
the same solution as the above equation, but do it in slightly different ways (see  docs).

#tip([Boundary data in LehrFEM])[
  `FixFlaggedSolutionComponents` (in both versions) requires a lambda function `std::pair<bool, double> selector(unsigned int dof_idx)`. It returns wether the dof is to be fixed and if so, the value it should be fixed to. 

  We get the boundary flags with

  ```cpp
  auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2);
  ```
  The second argument is the codim of entities we are interested in. To fix all dofs on the boundary (including those associated with edges and cells), give no second argument.

  *1. One function $g$ for the whole boundary*

  To get the function values of $g$ at the boundary nodes, wrap it in a `MeshFunction`:

  ```cpp
  auto mf_g = lf::mesh::utils::MeshFunctionGlobal(g);
  auto boundary_val = 
    lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g);
  ```
  This gives a `std::vector<std::pair<bool, scalar_t>>`. To get our selector:

  ```cpp
  auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
    return boundary_val[dof_idx];
  };
  ```
  If $g$ has different definitions on different parts of the boundary (e.g., $g=0$ on $Gamma_0$ and $g=1$ on $Gamma_1$), try to express $g$ as a lambda function with an if-else statement.

  *2. Constant value*

  The selector becomes

  ```cpp
  auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
    if (bd_flags[dof_idx]) {
      return std::make_pair(true, boundary_value);
    } else {
      return std::make_pair(false, 0.0); // value irrelevant
    }
  };
  ```
  *3. BC only on part of the boundary*

  Wee need to create our own `bd_flags`. Initialize it with default value `false`:
  
  ```cpp
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(mesh_p, false);
  ```
  Then loop over nodes and edges _separately_ and set `bd_flags` to `true` for the nodes/edges where the BC should be applied.
  ```cpp
  for (const auto& edge : fe_space->Mesh()->Entities(1)) {
    if (...) bd_flags(*edge) = true;
  }
  for (const auto& node : fe_space->Mesh()->Entities(2)) {...}
  ```
]
#pagebreak()

== Parametric Finite Element Methods
<sub:parametric-fem>
#lemma(number: "2.7.5.14", "Affine transformation of triangles")[
  For any triangle $K$ s.t. $lr(|K|) > 0$, there is a unique affine transformation
  $Phi_K (bxhat) = F_K bxhat + tau_K$, with $K = Phi_K (Khat)$ and
  $Khat$ the unit triangle. 
]<thm:affine-triangle-transformation>

This is nice since it allows us to transform integrals on an arbitrary triangle to ones
on the unit triangle and perform easy integration there. Additionally, $Phi_K$ can be
computed straightforwardly:
$ "Let " K "be a triangle with vertices " a^1 , a^2 , a^3. "Then" Phi_K (bxhat) = mat(delim: "[", a_1^2 - a_1^1, a_1^3 - a_1^2; a_2^2 - a_2^1, a_2^3 - a_2^2) bxhat + mat(delim: "[", a_1^1; a_2^1) $
#definition(number: "2.8.1.2",  "Pullback")[
  Given domains $Omega , mhat(Omega) subset bb(R)^d$ and a bijective mapping
  $Phi : mhat(Omega) arrow.r Omega$, the #emph[pullback] $Phi^(\*) u$ of a function
  $u : Omega arrow.r bb(R)$ is a function on $mhat(Omega)$ defined by
  $(Phi^(\*) u) (bxhat) := u (Phi (bxhat)) , bxhat in mhat(Omega)$.
]

An example for this is the affine transformation of triangles (@thm:affine-triangle-transformation) above. Note that in the following we will use $hat(x)$ for an
element that \"lives\" in a reference triangle (or quadrilateral). The
reference triangle has corners ${ (0 , 0) , (1 , 0) , (0 , 1) }$; the
reference quad is the unit square.

When dealing with #emph[parametric] finite element methods, we know that 
$ hat(b)_Khat^i = Phi_k ^(\*) b_K^i, $ where $b_K^i:K -> RR$ are the local basis functions on the concrete element $K$ and $hat(b)_Khat^i$ are the local basis functions on the reference element $Khat$ that matches $K$.

To make life easier, we call everything on the reference element *local* and everything on the concrete element *global*.

Note that all bilinear forms and linear forms in this course consist of
integrals, so we will sooner or later need to use quadrature to
approximate the integrals. But in the literature, quadrature rules are
defined on reference elements. Hence, we want a change of variables in
the integrals of the linear and bilinear forms such that we can apply
these quadrature rules. This change of variable is given by the pullback
function.

For the simple example of a mass matrix we get (by some multidimensional
analysis)
#neq($ (A_K)_(i j) = integral_K b_K^j (x) b_K^i (x) dif x = integral_(Khat) (Phi_K^(\*) b_K^j) (hat(x)) (Phi_K^(\*) b_K^i) (hat(x)) med sqrt(det (D Phi_K^top (hat(x)) D Phi_K (hat(x)))) dif x , $) <eq:mass-matrix-pullback>

where we have simply applied the pullback to both functions, switched to integration on $Khat$, and added an integration element. Similarly, for the diffusion matrix we get

#neq($ integral_K gradsub(bx) &b_K^i (x) dot.op gradsub(bx) b_K^j (x) dif x \
&= integral_(Khat) Phi_K^(\*) (gradsub(bx) b_K^j) (hat(x)) thin dot thin Phi_K^(\*) (gradsub(bx) b_K^i) (hat(x)) med sqrt(det (D Phi_K^top (hat(x)) D Phi_K (hat(x)))) dif x . $) <eq:diffusion-matrix-pullback>

by pulling back the global gradients (denoted by $gradsub(x)$) of the global shape functions.

Why is the integration element $det (D Phi_K^top (hat(x)) D Phi_K (hat(x)))$ more complicated than what we know from @thm:transformation-rule-for-integration? It is needed when $Omega$ and $mhat(Omega)$ do not live in the same space. If, for
example, $Omega subset bb(R)^3$ describes a 2D plane "living" in 3D space and
$mhat(Omega) subset bb(R)^2$, we will have
$D Phi_K in bb(R)^(3 times 2)$ and hence $det (Phi_K)$ is not defined, so the full term is used.
In the simpler case of square $Phi$, this simplifes to $det (Phi_K)$. For both cases, the term is provided by
#weblink("https://craffael.github.io/lehrfempp/classlf_1_1geometry_1_1_geometry.html#a80112cf5cfa9314cb44e61756299607d")[`lf::geometry::IntegrationElement`];.

The next thing which needs clarification is
$Phi_K^(\*) (gradsub(bx) b_K^i)$, pullback of the gradient.
As the gradient $gradsub(bx) b_K^j (x)$ depends on the shape of
$K$, it would be complicated to compute this term directly.

Therefore we use

#lemma(number: "2.8.3.10", "Transformation formula for gradients")[
  For differentiable $u : K arrow.r bb(R)$ and any diffeomorphism
  $Phi : Khat arrow.r K$ we have
  $ (gradsub(bxhat) Phi^(\*) u) (bxhat) = (D Phi (bxhat))^top underbrace([(gradsub(bx) u) (Phi(hat(x)))], (Phi^(\*) gradsub(bx)) u) (bxhat) $
  In words: on the left side we have *local gradient of the pullback* of $u$, on the right side we have *pullback of the global gradient* of $u$.
]
This is exactly what we need since the pullback of the global gradient occurs in @eq:diffusion-matrix-pullback and the local gradient of local shape functions is easy to compute.

Now we need to get rid of $(D Phi (hat(x)))^top$ on the right. Let's define $J=D Phi (hat(x))$. Since J might not be square, we cannot simply invert it. But we can multiply by $J (J^top J)^(- 1)$ to get 
$ Phi_K^(\*) (gradsub(bxhat) b_K^i) = J (J^top J)^(- 1) (gradsub(bxhat) hat(b)_Khat^i) $

$hat(b) (hat(x))$ is a reference basis function on the reference shape --- we can compute its gradient easily by hand.

Note that in the script, the case of square $J$ is assumed, such that this simplifies to $J^(- top) (gradsub(bxhat) hat(b)_Khat^i)$. The long term only matters
when we have (as above) $Omega$ and $mhat(Omega)$ not living in the same
space. The same remedy applies here, you can in any case use
#weblink("https://craffael.github.io/lehrfempp/classlf_1_1geometry_1_1_geometry.html#a7cb2b572966d7492522acb1b127cbbd0")[`lf::geometry::JacobianverseGramian`]
which will return $J^(- top)$ or $J (J^top J)^(- 1)$, respectively.

#block[
#strong[Bilinear Transformation for Quadrilaterals] Let
${ bold(a^1 \, a^2 \, a^3 \, a^4) }$ be the ordered corners of a
quadrilateral. Then
$ Phi_K (hat(x)) & = (1 - hat(x)_1) (1 - hat(x)_2) bold(a^1) + hat(x)_1 (1 - hat(x)_2) bold(a^2) + (1 - hat(x)_1) hat(x)_2 bold(a^3) + (1 - hat(x)_1) hat(x)_2 bold(a^4)\
 & = mat(delim: "[", alpha_1 + beta_1 hat(x)_1 + gamma_1 hat(x)_2 + delta_1 hat(x)_1 hat(x)_2; alpha_2 + beta_2 hat(x)_1 + gamma_2 hat(x)_2 + delta_2 hat(x)_1 hat(x)_2) $
with
$ mat(delim: "[", alpha_1; alpha_2) = bold(a^1) , mat(delim: "[", beta_1; beta_2) = bold(a^2 - a^1) , mat(delim: "[", gamma_1; gamma_2) = bold(a^4 - a^1) , mat(delim: "[", delta_1; delta_2) = bold(a^4 - a^3 - a^2 + a^1) $

]