# Function to plot the histogram of the posteriors on the final spin
# NKJ-M, 03.2016, based on earlier code

def af_posterior(fitsel, m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12, tag, outdir):

  import numpy as np
  import aux
  import nrutils_new as nr

  N_bin = 100 

  # Select the appropriate fit

  def fitaf(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, chi_pvec, phi12vec):
    if fitsel == 'HLZ':
      return nr.bbh_final_spin_non_precessing_Healyetal(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec))
    if fitsel == 'HLZinplane':
      return nr.bbh_final_spin_precessing_Healyetal_extension_Mf(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'HLZinplaneIMRPhenomP':
      return nr.bbh_final_spin_precessing_Healyetal_extension_Minit(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'Husaetal':
      return nr.bbh_final_spin_non_precessing_Husaetal(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec))
    if fitsel == 'IMRPhenomPv2':
      return nr.bbh_final_spin_precessing_Husaetal_extension_Minit_chi_p(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec), chi_pvec)
    if fitsel == 'IMRPhenomPtilts':
      return nr.bbh_final_spin_precessing_Husaetal_extension_Minit(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'BR':
      return nr.bbh_final_spin_precessing_Barausse_and_Rezzolla(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'BRaligned':
      return nr.bbh_final_spin_precessing_Barausse_and_Rezzolla(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec), 0.,0.,0.)

  tag = '%s_%s'%(tag,fitsel)

  # compute the posterior samples of the final mass and spin 
  chif = fitaf(m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12)

  print("Calculated posterior samples for chif\n")

  # Compute posterior distribution for chif

  P_chif, chif_bins = np.histogram(chif, bins=N_bin, normed=True)
  P_chif = P_chif.T
  print("Calculated histogram for chif\n")

  # Outputting the median and maP values and 68%, 90%, and 95% confidence intervals
                                                                
  chif_median = np.median(chif)
  chif_maP = aux.cen_from_edg(chif_bins)[np.nanargmax(P_chif)]

  chif_conf = aux.confidence_1d(P_chif, chif_bins)
                                                                                                                                                  
  chif_68min, chif_68max = chif_conf.edges(0.68)
  chif_90min, chif_90max = chif_conf.edges(0.90)
  chif_95min, chif_95max = chif_conf.edges(0.95)

  print("chif median: %.3f"%(chif_median))
  print("chif maP: %.3f"%(chif_maP))
  #print("chif 0.50 region: %.3f, %.3f"%(chif_50min, chif_50max))                                                                                     
                                                                                                                           
  print("chif 0.68 region: %.3f, %.3f"%(chif_68min, chif_68max))
  print("chif 0.90 region: %.3f, %.3f"%(chif_90min, chif_90max))
  print("chif 0.95 region: %.3f, %.3f\n"%(chif_95min, chif_95max))

  print("chif median \pm 0.90 edges: %.3f + %.3f - %.3f\n"%(chif_median, chif_90max - chif_median, chif_median - chif_90min))

  return chif
