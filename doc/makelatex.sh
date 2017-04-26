#!/bin/bash
cp ../README.rst .
for file in `ls *.rst`; do pre=${file%.rst}; (sed -e 's/<\([[:alpha:]]*\)>/<\1\.html>/g' $file | sed -e 's/https:\/\/bytebucket.org\/charade\/svengine\/wiki\///g' | rst2latex.py -) >$pre.tex; done;
for file in `ls *.tex`; do pre=${file%.tex}; pdflatex $file; done;
