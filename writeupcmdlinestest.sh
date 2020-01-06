#!/bin/bash

set -euo pipefail

cd $(dirname $0)

#cmd="sbatch --mem 4G --time 2:0:0 slurm.sh"; moreargs=""
#export SLURM_JOBID=123456; cmd="python"; moreargs=""
cmd="python"; moreargs=""
#cmd="srun --mem 4G --time 3:0:0 slurm.sh"; moreargs=""
#cmd="sbatch --mem 12G slurm.sh"; moreargs="onlyworkspace=1"

for scanrange in 101,-1,1:101,-0.02,0.02; do
  $cmd ./step9_runcombine.py fa3fa2fL1_EFT writeuptest productions=GEN_190908 scanranges=$scanrange scanfai=fa3 floatothers=0 runobs=0 usegs=0 $moreargs usesystematics=0
  $cmd ./step9_runcombine.py fa3fa2fL1_EFT writeuptest productions=GEN_190908 scanranges=$scanrange scanfai=fa2 floatothers=0 runobs=0 usegs=0 $moreargs usesystematics=0
  $cmd ./step9_runcombine.py fa3fa2fL1_EFT writeuptest productions=GEN_190908 scanranges=$scanrange scanfai=fL1 floatothers=0 runobs=0 usegs=0 $moreargs usesystematics=0
done
