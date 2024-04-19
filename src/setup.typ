#import "@preview/physica:0.9.3": eval, Order, Set
#import "boxes.typ": *

#let neq(content) = math.equation(
    block: true,
    numbering: "(1)",
    content,
)

#let colMath(x, color) = text(fill: color)[$#x$]
#let accentcolor = rgb(255, 0, 255)

#let weblink = (url, content) => {
  text(fill: blue)[#underline[#link(url)[#content]]]
}

#let unimportant = (
  fill: rgb("f6f6f6"), stroke: gray + 1.5pt
)

// Bold math for vectors
#let bx = $ bold(x) $
#let bu = $ bold(u) $
#let bw = $ bold(w) $
#let bA = $ bold(A) $
#let balpha = $ bold(alpha) $
#let bmu = $ bold(mu) $

// Wider hat
#let mhat = content => $ hat(content, size: #140%)$
#let Khat = $ hat(K, size: #120%)$
#let bxhat = $ hat(bx) $

#let grad = $bold("grad")thin $
#let gradsub = content => $bold("grad"_(#content)) $
#let div = $"div"thin $
#let openint(a,b) = $lr(\] #a, #b \[)$
#let fvH = $bold(cal(H))$

// #let eps = $ #h(0cm)text(epsilon, font: "Asana Math")med$
#let eps = $epsilon.alt$

#let this-template(doc) = [
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
  #set page(numbering: "1")
  #doc
]