#!/bin/bash
#cp ../README.rst HOME.rst
for file in `ls *.rst`; do pre=${file%.rst}; (sed -e 's/<\([[:alpha:]]*\)>/<\1\.html>/g' $file | rst2html.py -) >$pre.html; done;
ls *.html | zip SVEngine-doc.zip -@
