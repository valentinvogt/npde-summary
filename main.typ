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
      it.at("target"), numbering(el.numbering, ..counter(eq).at(el.location())),

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
#let date = datetime.today()

Last updated on #date.display("[month repr:long] [day], [year]")
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

#tip("Tips and practical advice")[
  from the TAs are highlighted like this.
]

= Basics
<ch:preface>
#v(-0.1cm)

#theorem(number: "0.3.1.19", "Cauchy-Schwarz Inequality")[
  If $a$ is a symmetric positive semi-definite bilinear form, then
#neq($ lr(|a (u , v)|) lt.eq a (u , u)^(1 / 2) a (v , v)^(1 / 2) $)
] <thm:cauchy-schwarz>

#v(-0.5cm)

#equation(number: "1.3.4.15", "Cauchy-Schwarz for Integrals")[
  #neq($ lr(|integral_Omega u(bx) v(bx) dif bx|) <= (integral_Omega |u(bx)|^2 dif bx)^(1 / 2) (integral_Omega |v(bx)|^2 dif bx)^(1 / 2) = norm(u)_(L^2 (Omega)) norm(v)_(L^2 (Omega)) $)
] <eq:cauchy-schwarz-integrals>

#v(-0.5cm)

#mybox("Norms")[
- $ bold("Supremum norm: ") norm(bold(u))_oo = norm(bold(u))_(L^oo (Omega)) := sup_(bx in Omega) " "norm(bold(u (x))) $
- $bold(L^2) bold("norm: ") norm(bold(u))_2 = norm(bold(u))_(L^2 (Omega)) := (integral_Omega norm(bold(u (x)))^2 dif x)^(1 / 2) $
]

#v(-0.5cm)

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