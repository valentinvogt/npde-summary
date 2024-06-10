#import "theorems.typ": *
#import "colors.typ": *

#let theorem = thmbox(
  "theorem", "Theorem", fill: brown-box,
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
  "lemma", "Lemma", fill: brown-box,
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
  "definition", "Definition", fill: brown-box,
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
  "equation", "Equation", fill: brown-box,
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
  "", "", fill: brown-box,
  namefmt: name => [
    #text(weight: "bold")[
      #h(-3pt) #name
    ]
  ],
  separator: linebreak()
).with(numbering: none)

// Modified from colorful-boxes:1.2.0
#let colorbox(title: none, color: none, radius: 2pt, width: auto, body) = {

  let strokeColor = tip-stroke
  let backgroundColor = page-color

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

#let tip-box = (title, content) => {
  colorbox(
    title: title,
    radius: 2pt,
  )[
    #v(-0.2cm)
    #content
    ]
}

#let subtle-box = (content, width: 100%) => {
  box(radius: 0.5cm, stroke: 1pt + text-color, inset: 0.5cm, width: width)[
    #content
  ]
}
