I hope everything here should be usable by everyone with no bugs, at least for the data and MC from 160720.

Step 2 can only be run on lxplus, unless you want to edit `helperstuff/samples.py` and change directory names to point to the CJLST trees.  The rest of the steps can probably be run somewhere else if you copy the step3 directory there.

To checkout:

`git clone --recursive git@github.com:hroskes/anomalouscouplings`

recursive is necessary, if you didn't do it, then after it's checked out run:

`git submodule update --init --recursive`


Then edit `helperstuff/config.py` with the directory where the repository is stored and the directory to store plots.

Then run:

```
python step1.py
python step2_adddiscriminants.py  #I generally run in parallel on 10 screens split across 2 lxplus machines
python step4_makejson.py  #quick, no need to run in parallel
python step6_maketemplates.py  #parallel again
python step8_combinesystematics.py   #quick
python step9_runcombine.py analysis foldername
    analysis: one of fa3 fa2 fL1
    foldername: gets appended to the analysis to form the folder name where the datacards are stored (-a for make_prop) and the folder where the plots are stored in plotsdir/limits (plotsdir specified in config)
    optional arguments, specified by name=value:
        plotname: default=limit, saves it as plotname.{png,eps,root,pdf}
        expectvalues: expectation scans for these values of fai, separated by commas, default is just 0
        legendposition: 4 numbers separated by commas to pass to the constructor of the TLegend.  Default: .2,.7,.6,.9
        channels: which channels to use for the scan, default is 2e2mu,4e,4mu (all of them)
        CLtextposition: x position along the 68% and 95% CL lines to put their labels.  Either left, right, or a float.  Default: left (equivalent to -1)
```
Running step9 will only run scans that haven't been run yet with the same analysis and foldername, so running after adjusting just the legend or CL text position is really quick.  The command used also gets saved as a .txt file with the same name as the plot so that it can be easily adjusted.  To remake all the plots using the same options after adjusting the plotting script, you can run
```
python step10_plottingutilities/replotlimits.py
```

To make plots, use all the python scripts in step10_plottingutilities.  Most of them take no arguments.  `projections.py` takes a while, and `makecontrolplots.py` takes a really long time.  `niceplots.py` is quick, but can only be run after `projections.py`.  `printrates.py` is pretty fast, it just prints the yields and number of observed events.

`limits.py` prints the 68% and 95% CL limits.  It has to be given the same arguments as `step9_runcombine.py`.

```
python limits.py analysis foldername
    optional arguments:
        plotname: has to be whatever was given to step9_runcombine.py.  Again limit is the default.
        format: latex or ppt, latex can be pasted into latex, ppt works with the microsoft equation editor
```
