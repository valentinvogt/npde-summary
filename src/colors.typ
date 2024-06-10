#let darkmode = "false"

#let page-color = white
#let text-color = black
#let brown-box = rgb(255, 236, 217)
#let tip-stroke = rgb(170, 170, 200)
#let unimportant-color = rgb("f6f6f6")
#let link-color = blue.darken(30%)

#if darkmode == "true" {
    page-color = rgb(39, 43, 50)
    text-color = white
    brown-box = rgb(80, 80, 80)
    tip-stroke = rgb(70, 70, 100)
    unimportant-color = gray
    link-color = blue.lighten(40%)
}