#!/usr/bin/env bash
if [[ -z $1 ]]; then
  echo "usage: runLongRanger lr.cfg"
  exit
fi

source $1
source $longranger_path
nohup nice longranger run --id=$id --sex=$sex --fastqs=$fastq_path --reference=$ref_path --indices=$indices --localcores=$localcores 1>$stdout 2>$stderr &
