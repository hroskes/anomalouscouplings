#!/bin/bash

set -euo pipefail

cd $(dirname $0)

cmd="sbatch slurm.sh"
#cmd=python

#found from decay only
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 floatothers=0
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 faiorder=fa3,fL1,fL1Zg,fa1,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 faiorder=fa3,fa1,fa2,fL1,fL1Zg
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 faiorder=fa3,fa2,fL1Zg,fa1,fL1

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 floatothers=0
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fL1,fa1,fL1Zg,fa3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fL1,fa3,fa1,fL1Zg
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fL1Zg,fa3,fa1,fL1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa1,fL1,fa3,fL1Zg
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa1,fL1Zg,fL1,fa3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fL1,fa1,fL1Zg
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fL1Zg,fa1,fL1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fa1,fL1Zg,fL1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:0.3,CMS_zz4l_fai4_relative:-0.8
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.1,CMS_zz4l_fai3_relative:0.1,CMS_zz4l_fai4_relative:-0.1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.1,CMS_zz4l_fai3_relative:0.2,CMS_zz4l_fai4_relative:-1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.1,CMS_zz4l_fai3_relative:0.3,CMS_zz4l_fai4_relative:-0.7
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.4,CMS_zz4l_fai3_relative:0.1,CMS_zz4l_fai4_relative:-1

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 floatothers=0
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai4_relative:-1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0.3,CMS_zz4l_fai2_relative:0.1,CMS_zz4l_fai4_relative:-1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0.4,CMS_zz4l_fai2_relative:-0.1,CMS_zz4l_fai4_relative:0.1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0.4,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai4_relative:-0.5
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa1,fa3,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa3,fa1,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa1,fL1Zg,fa2,fa3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa1,fL1Zg,fa3,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa1,fa3,fL1Zg,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa2,fa1,fL1Zg,fa3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa2,fa1,fa3,fL1Zg
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa2,fa3,fa1,fL1Zg
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa3,fa1,fL1Zg,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa3,fa2,fa1,fL1Zg

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg floatothers=0
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.1,CMS_zz4l_fai3_relative:0.3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai3_relative:0.2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai3_relative:0.4
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai2_relative:-1,CMS_zz4l_fai3_relative:0.4
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fL1,fa1,fa3,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fL1,fa3,fa1,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa1,fa2,fL1,fa3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa1,fa3,fL1,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa2,fa1,fa3,fL1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa2,fa3,fa1,fL1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa3,fL1,fa1,fa2
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa3,fa2,fa1,fL1

#more for decay

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa1,fa2,fa3
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa3,fa2,fa1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa2,fa3,fa1
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa2,fa1,fa3

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:-0.025,CMS_zz4l_fai4_relative:0.9

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg yieldsystematics scanranges=101,-1,1:101,-0.02,0.02 plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fL1,fa1,fL1Zg setparametersforgrid=CMS_zz4l_fa1_relative:0.18,CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:0
