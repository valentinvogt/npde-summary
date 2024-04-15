#import "theorems.typ": *

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

// Modified from colorful-boxes:1.2.0
#let colorbox(title: none, color: none, radius: 2pt, width: auto, body) = {

  let strokeColor = luma(70)
  let backgroundColor = white

  if color == "red" {
    strokeColor = rgb(237, 32, 84)
    backgroundColor = rgb(253, 228, 224)
  } else if color == "green" {
    strokeColor = rgb(102, 174, 62)
    backgroundColor = rgb(235, 244, 222)
  } else if color == "blue" {
    strokeColor = rgb(29, 144, 208)
    backgroundColor = rgb(232, 246, 253)
  }

  return box(
    fill: backgroundColor,
    stroke: 2pt + strokeColor,
    radius: radius,
    width: width
  )[
    #if (title != "") {
      block(
        fill: strokeColor, 
        inset: 8pt,
        radius: (top-left: radius, bottom-right: radius),
      )[
        #text(fill: white, weight: "bold")[#title]
      ]
    }
    #block(
      width: 100%,
      inset: if title != ""{
        (x: 8pt, bottom: 8pt)
      } else {
        (x: 8pt, bottom: 8pt, top: 13pt)
      }
    )[
      #body
    ]
  ]
}

#let tip_box = (title, content) => {
  colorbox(
    title: title,
    color: "blue",
    radius: 2pt,
  )[
    #v(-0.2cm)
    #content
    ]
}

#let subtle-box = (content, width: 100%) => {
  box(radius: 0.5cm, stroke: 1pt, inset: 0.5cm, width: width)[
    #content
  ]
}