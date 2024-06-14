#!/bin/bash

echo "converting images"
for i in *.pdf; do
#   pdftoppm -png -rx 300 -ry 500 $i > ${i%.pdf*}.png
   pdftoppm -png $i > ${i%.pdf*}.png
done

