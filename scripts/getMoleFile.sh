#!/bin/bash
#to use this script need to set $cfgfile and $moleFile
source $cfgfile
log=$moleFile.log.getMoleFile.$mydate
maxSize=$(cat $moleFile | wc -l)
maxSize=$(expr $(($maxPlot>$maxSize?$maxSize:$maxPlot)) - 1)
for i in `seq 1 $maxSize`; do
  getMolCMD="$ZOOMXPATH/../scripts/getMolecule.sh $moleFile $outputPf/$caseBam $i $moleFile.Plot >>$log 2>&1"
  echo $getMolCMD | bash
  #plotMolCMD="$ZOOMXPATH/../scripts/plotMolecule.R $isUCSC $moleFile.Plot_$i >>$log 2>&1"
  #echo $plotMolCMD | bash
done
