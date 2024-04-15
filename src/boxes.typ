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
