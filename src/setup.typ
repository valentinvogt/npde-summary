#import "@preview/physica:0.9.3": eval, Order
#import "boxes.typ": *

#let neq(content) = math.equation(
    block: true,
    numbering: "(1)",
    content,
)

#let colMath(x, color) = text(fill: color)[$#x$]
#let my_accent = rgb(255, 0, 255)

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

#let mhat = content => $ hat(content, size: #140%)$
#let Khat = $ hat(K, size: #120%)$
#let bxhat = $ hat(bx) $

#let grad = $bold("grad")thin $
#let gradsub = content => $bold("grad"_(#content)) $
#let div = $"div"thin $
