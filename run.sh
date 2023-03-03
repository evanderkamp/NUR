#!/bin/bash

echo "Run handin"

echo "Download txt file to use in script"
if [ ! -e Vandermonde.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

# First exercise
echo "Run the first script ..."
python3 NUR_handin1.py

# Second exercise
echo "Run the second script ..."
python3 NUR_handin1Q2.py


echo "Generating the pdf"

pdflatex Handin1.tex
bibtex Handin1.aux
pdflatex Handin1.tex
pdflatex Handin1.tex
