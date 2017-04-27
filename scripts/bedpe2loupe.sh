#!/usr/bin/env bash
if [[ -z $1 ]]; then
  echo "usage: bedpe2loupe.sh moleFile [chr]"
  echo "example: bedpe2loupe.sh NA12878.wg12.chr22.bam.dis.E40.refined" 
  echo "  convert bedpe coordinates to loupe coordinates. add [chr] to prefix if specified."
  exit
fi 
if [[ -z $2 ]]; then
  awk '{print $1":"$2"-"$3";"$4":"$5"-"$6"\t"$7"\t"$8";"$9";"$10}' $1
else
  awk '{print "chr"$1":"$2"-"$3";chr"$4":"$5"-"$6"\t"$7"\t"$8";"$9";"$10}' $1
fi
