#!/usr/bin/env bash
rm -f testDependency.log
echo "testing R dependencies ..."
#r_deps=(ak47 optparse rjson RColorBrewer IRanges GenomicRanges ggbio GenomeInfoDb)
r_deps=(optparse rjson RColorBrewer IRanges GenomicRanges ggbio GenomeInfoDb)
failed_r_deps="Failed R Dependencies: "
for d in ${r_deps[*]}; do 
  echo -n $d:
  cmd="echo \"library($d)\" | R --vanilla 1>>testDependency.log 2>&1"
  echo $cmd | bash
  if [[ ! $? -eq 0 ]]; then
    echo [FAILED]
    echo Failed Test: $cmd
    failed_r_deps="$failed_r_deps $d"
  else
    echo [OK]
  fi
done
echo "$failed_r_deps"
echo "testing Python dependencies ..."
#python_deps=(ak47 pysam pybedtools pandas tables scipy)
python_deps=(pysam pybedtools pandas tables scipy)
failed_python_deps="Failed Python Dependencies: "
for d in ${python_deps[*]}; do 
  echo -n $d:
  cmd="python -c \"import $d\" 1>>testDependency.log 2>&1"
  echo $cmd | bash
  if [[ ! $? -eq 0 ]]; then
    echo [FAILED]
    echo Failed Test: $cmd
    failed_python_deps="$failed_python_deps $d"
  else
    echo [OK]
  fi
done
echo "$failed_python_deps"
