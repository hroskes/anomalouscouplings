#!/bin/bash

set -euo pipefail

cd $(dirname $0)

cmd="sbatch --mem 12G slurm.sh"; moreargs=
#cmd=python; moreargs=
#cmd="sbatch --mem 12G slurm.sh"; moreargs="onlyworkspace=1"

for scanrange in 26,-1,1 26,-0.02,0.02; do
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 floatothers=0 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 faiorder=fa3,fL1,fL1Zg,fa1,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 faiorder=fa3,fa1,fa2,fL1,fL1Zg $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa3 faiorder=fa3,fa2,fL1Zg,fa1,fL1 $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 floatothers=0 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fL1,fa1,fL1Zg,fa3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fL1,fa3,fa1,fL1Zg $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fL1Zg,fa3,fa1,fL1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa1,fL1,fa3,fL1Zg $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa1,fL1Zg,fL1,fa3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fL1,fa1,fL1Zg $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fL1Zg,fa1,fL1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fa1,fL1Zg,fL1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:0.3,CMS_zz4l_fai4_relative:-0.8 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.1,CMS_zz4l_fai3_relative:0.1,CMS_zz4l_fai4_relative:-0.1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.1,CMS_zz4l_fai3_relative:0.2,CMS_zz4l_fai4_relative:-1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.1,CMS_zz4l_fai3_relative:0.3,CMS_zz4l_fai4_relative:-0.7 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.4,CMS_zz4l_fai3_relative:0.1,CMS_zz4l_fai4_relative:-1 $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 floatothers=0 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai4_relative:-1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0.3,CMS_zz4l_fai2_relative:0.1,CMS_zz4l_fai4_relative:-1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0.4,CMS_zz4l_fai2_relative:-0.1,CMS_zz4l_fai4_relative:0.1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0.4,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai4_relative:-0.5 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa1,fa3,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa3,fa1,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa1,fL1Zg,fa2,fa3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa1,fL1Zg,fa3,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa1,fa3,fL1Zg,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa2,fa1,fL1Zg,fa3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa2,fa1,fa3,fL1Zg $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa2,fa3,fa1,fL1Zg $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa3,fa1,fL1Zg,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fa3,fa2,fa1,fL1Zg $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg floatothers=0 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.1,CMS_zz4l_fai3_relative:0.3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai3_relative:0.2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.2,CMS_zz4l_fai3_relative:0.4 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai2_relative:-1,CMS_zz4l_fai3_relative:0.4 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fL1,fa1,fa3,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fL1,fa3,fa1,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa1,fa2,fL1,fa3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa1,fa3,fL1,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa2,fa1,fa3,fL1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa2,fa3,fa1,fL1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa3,fL1,fa1,fa2 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa3,fa2,fa1,fL1 $moreargs

#more for decay

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa1,fa2,fa3 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa3,fa2,fa1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa2,fa3,fa1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 faiorder=fL1,fL1Zg,fa2,fa1,fa3 $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:-0.025,CMS_zz4l_fai4_relative:0.9 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 faiorder=fa2,fa3,fL1,fa1,fL1Zg setparametersforgrid=CMS_zz4l_fa1_relative:0.18,CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:0 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0.1,CMS_zz4l_fai3_relative:0.3,CMS_zz4l_fai4_relative:-1 $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fL1,fa1,fa3,fa2 setparametersforgrid=CMS_zz4l_fai2_relative:-1,CMS_zz4l_fai3_relative:0.4 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg faiorder=fL1Zg,fa3,fa2,fa1,fL1 setparametersforgrid=CMS_zz4l_fai2_relative:-0.6,CMS_zz4l_fa1_relative:0 $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:0.2,CMS_zz4l_fai4_relative:-1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:0.23,CMS_zz4l_fai4_relative:-1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:0.26,CMS_zz4l_fai4_relative:-1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:-0.01,CMS_zz4l_fai4_relative:0.6 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai3,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fa2 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai3_relative:-0.016,CMS_zz4l_fai4_relative:0.74 $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:0.006,CMS_zz4l_fai4_relative:0.89 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:0.005,CMS_zz4l_fai4_relative:0.82 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:0.005,CMS_zz4l_fai4_relative:0.80 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:0.004,CMS_zz4l_fai4_relative:0.74 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:0.003,CMS_zz4l_fai4_relative:0.68 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai4,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1 setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:0.002,CMS_zz4l_fai4_relative:0.62 $moreargs

$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0.11,CMS_zz4l_fai2_relative:-0.03,CMS_zz4l_fai3_relative:1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0.12,CMS_zz4l_fai2_relative:-0.02,CMS_zz4l_fai3_relative:1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.3,CMS_zz4l_fai3_relative:1 $moreargs
$cmd ./step9_runcombine.py fa3fa2fL1fL1Zg_morecategories newyields scanranges=$scanrange plotnuisances=CMS_zz4l_fai1,CMS_zz4l_fai2,CMS_zz4l_fai3,CMS_zz4l_fa1 useNLLandNLL0=0 scanfai=fL1Zg setparametersforgrid=CMS_zz4l_fai1_relative:0,CMS_zz4l_fai2_relative:-0.24,CMS_zz4l_fai3_relative:1 $moreargs
done

