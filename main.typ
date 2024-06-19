#import "src/theorems.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "src/setup.typ": *
#show: this-template

#align(
  center,
)[
  = Numerical Methods for PDEs --- TA Summary

  #v(0.5cm)

  2023 version created by
  #link("mailto:jbachmann@student.ethz.ch")[Jonas Bachmann],
  #link("mailto:pfischill@student.ethz.ch")[Paul Fischill],
  #link("mailto:samrusso@student.ethz.ch")[Samuel Russo] and
  #link("mailto:grafn@student.ethz.ch")[Nico Graf]. \

  2024 version updated by
  #link("mailto:vogtva@student.ethz.ch")[Valentin Vogt],
  #link("mailto:luwirth@student.ethz.ch")[Luis Wirth].

  #v(0.5cm)
  #let date = datetime.today()

  Last updated on #date.display("[year]-[month]-[day]").
]

= About
#v(-0.1cm)

#mybox("Theorems, definitions and equations")[
  from the lecture notes come in boxes like this one.
]
#v(-0.5cm)

#mybox("Less important results", ..unimportant)[
  that are given for context or completeness look like this.
]
#v(-0.3cm)

#tip-box("Tips and practical advice")[
  from the TAs are highlighted like this.
]

= Basics
#v(-0.1cm)

#theorem(number: "0.3.1.19", [Cauchy--Schwarz Inequality])[
  If $a$ is a symmetric positive semi-definite bilinear form, then
  #neq($ lr(|a (u , v)|) lt.eq a (u , u)^(1 / 2) a (v , v)^(1 / 2) $)
] <thm:cauchy-schwarz>

#equation(
  number: "1.3.4.15", [Cauchy--Schwarz for Integrals],
)[
  #neq(
    $ lr(|integral_Omega u(bx) v(bx) dif bx|) <= (integral_Omega |u(bx)|^2 dif bx)^(1 / 2) (integral_Omega |v(bx)|^2 dif bx)^(1 / 2) = norm(u)_(L^2 (Omega)) norm(v)_(L^2 (Omega)) $,
  )
] <eq:cauchy-schwarz-integrals>

#mybox(
  "Norms",
)[
  - $ bold("Supremum norm: ") norm(bold(u))_oo = norm(bold(u))_(L^oo (Omega)) := sup_(bx in Omega) " "norm(bold(u (x))) $
  - $bold(L^2) bold("norm: ") norm(bold(u))_2 = norm(bold(u))_(L^2 (Omega)) := (integral_Omega norm(bold(u (x)))^2 dif x)^(1 / 2) $
]

#mybox(
  [Barycentric coordinate functions],
)[
  The barycentric coordinate functions $lambda_i$ on a triangle with vertices $ba_1, ba_2, ba_3$ are linear functions satisfying the *cardinal property*
  #neq(
    $ lambda_i (ba_j) = delta_(i j) = cases(
      1 & quad upright("if") i = j,
      0 & quad upright("else")
    ) $
  )<eq:barycentric-cardinal-property>
  
  On the unit triangle, whose vertices are $ba_1 = (0, 0), ba_2 = (1, 0), ba_3 = (0, 1)$, the barycentric coordinate functions are 
  #neq(
    $ lambda_1 (bx) &= 1 - x_1 - x_2 \
      lambda_2 (bx) &= x_1 \
      lambda_3 (bx) &= x_2 $
  )
]
#theorem(
  number: "0.3.2.31", "Transformation rule for Integration",
)[
  Given two domains $Omega , mhat(Omega)$ and a continuous, differentiable mapping $Phi : mhat(Omega) arrow.r Omega$
  #neq(
    $ integral_Omega f (bx) dif bx = integral_(mhat(Omega)) f (Phi (hat(x))) lr(|det D Phi (hat(x))|) dif bold(hat(x)) $,
  )
] <thm:transf-rule-integration>

#set heading(numbering: "1.1")

#pagebreak()
#include "src/chapters/01.typ"

#pagebreak()
#include "src/chapters/02.typ"

#pagebreak()
#include "src/chapters/03.typ"

#pagebreak()
#counter(heading).update(4)
#include "src/chapters/05.typ"

#pagebreak()
#counter(heading).update(8)
#include "src/chapters/09.typ"

#pagebreak()
#include "src/chapters/10.typ"

#pagebreak()
#include "src/chapters/11.typ"

#pagebreak()
#include "src/chapters/12.typ"
