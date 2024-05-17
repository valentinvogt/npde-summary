#!/usr/bin/env sh

typst watch main.typ out/summary.pdf &

while :
do
  zathura out/summary.pdf
done
