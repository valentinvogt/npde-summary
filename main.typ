#import "src/theorems.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "src/setup.typ": *
#show: thmrules
#show ref: it => {
  let eq = math.equation
  let el = it.element
  if el != none and el.func() == eq {
    // Override equation references.
    link(
      it.at("target"),
      numbering(
      el.numbering,
      ..counter(eq).at(el.location())),

    )
[]
    
  } else {
    it
  }
}


#align(center)[
= Numerical Methods for PDEs — TA Summary

#v(0.5cm)

Created by #weblink("mailto:jbachmann@student.ethz.ch")[Jonas Bachmann];,
#weblink("mailto:pfischill@student.ethz.ch")[Paul Fischill];,
#weblink("mailto:samrusso@student.ethz.ch")[Samuel Russo] and
#weblink("mailto:grafn@student.ethz.ch")[Nico Graf] \

Ported to Typst by #weblink("mailto:vogtva@student.ethz.ch")[Valentin Vogt] \
#v(0.5cm)

Last updated on April 11, 2024.
]

= Basics
<ch:preface>

#theorem(number: "0.3.1.19", "Cauchy-Schwarz Inequality")[
  If $a$ is a symmetric
positive semi-definite bilinear form, then
#neq($ lr(|a (u , v)|) lt.eq a (u , u)^(1 / 2) a (v , v)^(1 / 2) $)
] <thm:cauchy-schwarz>
#mybox("Norms")[
- $ bold("Supremum norm: ") norm(bold(u))_oo = norm(bold(u))_(L^oo (Omega)) := sup_(bx in Omega) " "norm(bold(u (x))) $
- $bold(L^2) bold("norm: ") norm(bold(u))_2 = norm(bold(u))_(L^2 (Omega)) := (integral_Omega norm(bold(u (x)))^2 dif x)^(1 / 2) $
]

#theorem(number: "0.3.2.31", "Transformation rule for Integration")[
  Given two domains $Omega , mhat(Omega)$ and a continuous, differentiable mapping $bold(Phi) : mhat(Omega) arrow.r Omega$
  #neq($ integral_Omega f (bx) dif bx = integral_(mhat(Omega)) f bold((Phi (hat(x)))) lr(|det D bold(Phi (hat(x)))|) dif bold(hat(x)) $)
] <thm:transformation-rule-for-integration>

#set heading(numbering: "1.1")

#pagebreak()
#include "chapters/01.typ"
#pagebreak()
#include "chapters/02.typ"
#pagebreak()
#include "chapters/03.typ"
#pagebreak()
#counter(heading).update(4)
#include "chapters/05.typ"
#counter(heading).update(8)
#include "chapters/09.typ"
#pagebreak()
#include "chapters/10.typ"
#pagebreak()
#include "chapters/11.typ"