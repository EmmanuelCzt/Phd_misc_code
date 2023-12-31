#!/bin/sh

# Convert all arguments (assumed SVG) to a jpg
# Requires Inkscape and ImageMagick 6.8 (doesn't work with 6.6.9)
# Source : https://gist.github.com/matsen/4263955

for i in $@; do
  BN=$(basename $i .svg)
  inkscape --without-gui --export-png="$BN.png" --export-dpi 300 $i
  convert -compress LZW -alpha remove $BN.png $BN.jpg
  mogrify -alpha off $BN.jpg
  rm $BN.png
done
