#import "@preview/physica:0.9.3": eval, Order, Set
#import "boxes.typ": *
#import "@preview/natrix:0.1.0": *
#import "colors.typ": *

#let neq(content) = math.equation(
    block: true,
    numbering: "(1)",
    content,
)

#let colMath(x, color) = text(fill: color)[$#x$]
#let accentcolor = rgb(255, 0, 255)

#let unimportant = (
  fill: unimportant-color, stroke: unimportant-color.darken(20%) + 1.5pt,
)

// Bold math for vectors
#let ba = $ bold(a) $
#let bm = $ bold(m) $
#let bx = $ bold(x) $
#let bu = $ bold(u) $
#let bv = $ bold(v) $
#let bw = $ bold(w) $
#let bA = $ bold(A) $
#let balpha = $ bold(alpha) $
#let bmu = $ bold(mu) $
#let bH = $ bold(H) $
#let jac = $ upright("D") $
// Wider hat
#let mhat = content => $ hat(content, size: #140%)$
#let Khat = $ hat(K, size: #120%)$
#let bxhat = $ hat(bx) $
#let msh = $ cal(M) $

#let grad = $bold("grad")thin $
#let gradsub = content => $bold("grad"_(#content)) $
#let div = $"div"thin $
#let openint(a,b) = $lr(\] #a, #b \[)$
#let fvH = $bold(cal(H))$
#let curl = $bold("curl")thin $

#let recop = $upright("R")_M arrow(bmu)$
#let eps = $epsilon.alt$
#let argmin = math.op("arg min", limits: true)
#let argmax = math.op("arg max", limits: true)

#let this-template(doc) = [
  #show: thmrules
  #show link: set text(fill: link-color)

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

  #set text(text-color)
  #set page(fill: page-color)
  #doc
]
