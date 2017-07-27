#!/usr/bin/env bash

vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}

testvercomp () {
    vercomp $1 $2
    case $? in
        0) op='=';;
        1) op='>';;
        2) op='<';;
    esac
    if [[ $op != $3 ]]
    then
        echo "FAIL: Expected '$3', Actual '$op', Arg1 '$1', Arg2 '$2'"
    else
        echo "Pass: '$1 $op $2'"
    fi
}

rm -f test_dep.log
c_deps=(samtools bedtools)
failed_c_deps="Failed C Dependencies: "
echo "testing C dependencies ..."
for d in ${c_deps[*]}; do 
  echo -n $d:
  case $d in
    "samtools") ver=`samtools 2>&1 | grep -o "Version: [0-9.]\+" | grep -o "[0-9.]\+"`; 
    echo -n $ver; vercomp $ver "0.1.18"; op=$?;;
    "bwa") ver=`bwa 2>&1 | grep -o "Version: [0-9.]\+" | grep -o "[0-9.]\+"`; 
    echo -n $ver; vercomp $ver "0.7.0"; op=$?;;
    "bedtools") ver=`bedtools --version 2>&1 | grep -o "v[0-9.]\+" | grep -o "[0-9.]\+"`; 
    echo -n $ver; vercomp $ver "2.0.0"; op=$?;;
  esac
  if [[ ! $op -eq 1 ]]; then
    echo [FAILED]
    failed_c_deps="$failed_c_deps $d"
  else
    echo [OK]
  fi
done
echo "$failed_c_deps"

echo "testing R dependencies ..."
r_deps=(optparse rjson ggplot2 RColorBrewer cowplot igraph IRanges GenomicRanges ggbio GenomeInfoDb BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg38)
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
