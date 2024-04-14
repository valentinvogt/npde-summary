#import "theorems.typ": *
#import "@preview/colorful-boxes:1.2.0": *

#let theorem = thmbox(
  "theorem", "Theorem", fill: rgb("#ffecd9"),
  // bodyfmt: body => [
  //   #body 2
  // ]
  titlefmt: title => [
    #text(weight: "bold")[
      #title:
    ]
  ],
  namefmt: name => [
    #text(weight: "bold")[
      #name
    ]
  ],
  separator: linebreak(),
  supplement: "Theorem"
)

#let lemma = thmbox(
  "lemma", "Lemma", fill: rgb("#ffecd9"),
  titlefmt: title => [
    #text(weight: "bold")[
      #title:
    ]
  ],
  namefmt: name => [
    #text(weight: "bold")[
      #name
    ]
  ],
  separator: linebreak(),
  supplement: "Lemma"
)

#let definition = thmbox(
  "definition", "Definition", fill: rgb("#ffecd9"),
  titlefmt: title => [
    #text(weight: "bold")[
      #title:
    ]
  ],
  namefmt: name => [
    #text(weight: "bold")[
      #name
    ]
  ],
  separator: linebreak(),
  supplement: "Definition"
)

#let equation = thmbox(
  "equation", "Equation", fill: rgb("#ffecd9"),
  titlefmt: title => [
    #text(weight: "bold")[
      #title
    ]
  ],
  namefmt: name => [
    #text(weight: "bold")[
      (#name)
    ]
  ],
  separator: linebreak(),
  supplement: "Equation"
)

#let mybox = thmbox(
  "", "", fill: rgb("#ffecd9"),
  namefmt: name => [
    #text(weight: "bold")[
      #h(-3pt) #name
    ]
  ],
  separator: linebreak()
).with(numbering: none)

#let neq(content) = math.equation(
    block: true,
    numbering: "(1)",
    content,
)

#let colMath(x, color) = text(fill: color)[$#x$]
#let my_accent = rgb(255, 0, 255)

#show ref: it => {
  let eq = math.equation
  let el = it.element
  if el != none and el.func() == eq {
    // Override equation references.
    link(
      // it.at("target"),
      numbering(
      el.numbering,
      ..counter(eq).at(el.location()))
    )[]
    
  } else {
    it
  }
}

#let weblink = (url, content) => {
  text(fill: blue)[#underline[#link(url)[#content]]]
}

#let tip = (title, content) => {
  colorbox(
    title: title,
    color: "blue",
    radius: 2pt,
  )[
    #content
  ]
}

#let subtle-box = (content, width: 100%) => {
  box(radius: 0.5cm, stroke: 1pt, inset: 0.5cm, width: width)[
    #content
  ]
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