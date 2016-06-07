import ROOT

"""
Hi Heshy,
 
Let me include a portion of the email I sent to Aurelijus Rinkevicius some time ago at the end. Check Point e in particular. When you reweight from one hypothesis to another, you must make sure the sum of weights in 2l2l' for LHE-level events are the same; you cannot allow them to fluctuate, otherwise you will be dependent on the statistics of your sample more (the fact that these sums have to be te same is an additional constraint to improve your efficiency prediction). This number is essentially supposed to be the 2e2mu (*3 since we include taus) cross section*BR (JHUGen adjusted to YR4, with adjustment to 4.07 MeV width). The way we did was to normalize everything after reweighting to
a) The number of events in each sample
b) Renormalize (a) to have the same sum of weighted 2l2l' numbers for each particular hypothesis
c) Renormalize (b) such that 2l2l' (index[f=2] below) xsec of SM (index[p=0] below) is the same as YR3/4 prediction.
 
double total_weighted = NGenTotal_weighted[smp][p];
double total_unweighted = NGenTotal_unweighted[smp];
flavor_contribution[f] += NGenTotal_weighted_perFlavor[smp][p][f] * total_unweighted/total_weighted;
sum_nevents += NGenTotal_weighted_perFlavor[smp][p][f] * total_unweighted/total_weighted;
fGen_BRweight_perFlavor[p][f] = flavor_contribution[f] / sum_nevents;
fGen_BRweight[p] = fGen_BRweight_perFlavor[0][2]/fGen_BRweight_perFlavor[p][2];
MC_weight_spin0[p] *= ( NGenTotal_unweighted[smp]/NGenTotal_weighted[smp][p]*fGen_BRweight[p] );
float myxsec = MC_weight_xsec * MC_weight_spin0[p]  / sum_NGenTotal_unweighted ;
MC_CV_weight[p] = myxsec*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT;
 
After this procedure, all of your hypothesis will be normalized in 2l2l' channel to SM xsec, and you will just have to multiply by the mixture xsec/SM xsec from JHUGen (this I why I kept saying that having these numbers precisely is very important for mixture templates).

(...)

Ulascan
"""

def normalizeweights():
   samples = [Sample("ggH", hypothesis) for hypothesis in hypotheses]
   f = OrderedDict()
   t = OrderedDict()
   h = OrderedDict()
   newf = OrderedDict()
   newt = OrderedDict()
   for sample in samples:
       f[sample] = ROOT.TFile(sample.withdiscriminantsfile())
       t[sample] = f[sample].candTree
       h[sample] = f[sample].nevents

       newf[sample] = ROOT.TFile(sample.normalizedfile(), "recreate")
       newt[sample] = t[sample].Clone(0)
