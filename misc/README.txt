Using this branch of CombinedLimit:
https://github.com/hroskes/HiggsAnalysis-CombinedLimit/tree/HiggsJCPWidthLifetime
(or more precisely, this snapshot:
https://github.com/hroskes/HiggsAnalysis-CombinedLimit/tree/943fabe1dceae41a81c67b779f304540d36f11b1)

combineCards.py *.txt > hzz4l_4l_lumi35.8671.txt    #or *Untagged*.txt to scan only in the untagged category, etc.
text2workspace.py -m 125 hzz4l_4l_lumi35.8671.txt -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
                  --PO verbose --PO allowPMF -o workspace_lumi35.8671.root -v 7
combine -M MultiDimFit workspace_lumi35.8671.root --algo=grid --points 100 \
        --setPhysicsModelParameterRanges CMS_zz4l_fai1=-1.0,1.0 -m 125 -n $1_exp_0.0_lumi35.8671 \
        -t -1 --setPhysicsModelParameters r_ffH=1,r_VVH=1,CMS_zz4l_fai1=0.0 -V -v 3 --saveNLL \
        -S 1 --saveSpecifiedFunc=r_VVH,muV_scaled,r_ffH,muf_scaled --saveSpecifiedNuis all
