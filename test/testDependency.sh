#!/usr/bin/env bash
rm -f testDependency.log
echo "testing R dependencies ..."
r_deps=(optparse rjson RColorBrewer IRanges GenomicRanges ggbio GenomeInfoDb)
for d in ${r_deps[*]}; do 
  echo -n $d:
  cmd="echo \"library($d)\" | R --vanilla 1>>testDependency.log 2>&1"
  echo $cmd | bash
  if [[ ! $? -eq 0 ]]; then
    echo [FAILED]
    echo Failed Test: $cmd
    exit
  fi
  echo [OK]
done
echo "testing Python dependencies ..."
python_deps=(pysam pybedtools pandas tables scipy)
for d in ${python_deps[*]}; do 
  echo -n $d:
  cmd="python -c \"import $d\" 1>>testDependency.log 2>&1"
  echo $cmd | bash
  if [[ ! $? -eq 0 ]]; then
    echo [FAILED]
    echo Failed Test: $cmd
    exit
  fi
  echo [OK]
done
