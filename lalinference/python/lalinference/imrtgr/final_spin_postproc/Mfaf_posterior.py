# Function to plot the histogram of the posteriors on the final mass and spin
# NKJ-M, 03.2016, based on earlier code

def Mfaf_posterior(fitsel, m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12, tag, outdir):

  import numpy as np
  import aux
  import nr_fits_Mf_af as nrf

  N_bin = 100 

  # Select the appropriate fit

  def fitMfaf(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, chi_pvec, phi12vec):
    if fitsel == 'HLZ':
      return nrf.bbh_final_mass_and_spin_non_precessing(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec))
    if fitsel == 'HLZinplane':
      return nrf.bbh_final_mass_and_spin_HLZ_extension_precessing(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'HLZinplaneIMRPhenomP':
      return nrf.bbh_final_mass_and_spin_HLZ_extension_precessing_IMRPhenomP(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'Husaetal':
      return nrf.bbh_final_mass_and_spin_non_precessing_Husaetal(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec))
    if fitsel == 'IMRPhenomPv2':
      return nrf.bbh_final_mass_and_spin_precessing_IMRPhenomPv2(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec), chi_pvec)
    if fitsel == 'IMRPhenomPtilts':
      return nrf.bbh_final_mass_and_spin_precessing_IMRPhenomP_tilts(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'BR':
      return nrf.bbh_final_mass_and_spin_precessing_HLZ_mass_Barausse_and_Rezzolla_spin(m1vec, m2vec, chi1vec, chi2vec, tilt1vec, tilt2vec, phi12vec)
    if fitsel == 'BRaligned':
      return nrf.bbh_final_mass_and_spin_precessing_HLZ_mass_Barausse_and_Rezzolla_spin(m1vec, m2vec, chi1vec*np.cos(tilt1vec), chi2vec*np.cos(tilt2vec), 0.,0.,0.)

  tag = '%s_%s'%(tag,fitsel)

  # compute the posterior samples of the final mass and spin 
  Mf, chif = fitMfaf(m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12)

  print("Calculated posterior samples for Mf and chif\n")

  # Compute posterior distribution for Mf and chif

  P_Mf, Mf_bins = np.histogram(Mf, bins=N_bin, normed=True)
  P_Mf = P_Mf.T
  P_chif, chif_bins = np.histogram(chif, bins=N_bin, normed=True)
  P_chif = P_chif.T
  print("Calculated histogram for Mf and chif\n")

  # Outputting the median and maP values and 68%, 90%, and 95% confidence intervals for the original posterior                                                                       
  Mf_median = np.median(Mf)
  Mf_maP = aux.cen_from_edg(Mf_bins)[np.nanargmax(P_Mf)]

  Mf_conf = aux.confidence_1d(P_Mf, Mf_bins)

  #Mf_50min, Mf_50max = Mf_conf.edges(0.5)                                                                                                                                           
  Mf_68min, Mf_68max = Mf_conf.edges(0.68)
  Mf_90min, Mf_90max = Mf_conf.edges(0.90)
  Mf_95min, Mf_95max = Mf_conf.edges(0.95)
  
  print("Original values (all in solar masses):")
  print("Mf median: %.3f"%(Mf_median))
  print("Mf maP: %.3f"%(Mf_maP))
  #print("Mf 0.50 region: %.3f, %.3f"%(Mf_50min, Mf_50max))                                                                                                                        
  print("Mf 0.68 region: %.3f, %.3f"%(Mf_68min, Mf_68max))
  print("Mf 0.90 region: %.3f, %.3f"%(Mf_90min, Mf_90max))
  print("Mf 0.95 region: %.3f, %.3f\n"%(Mf_95min, Mf_95max))

  print("Mf median \pm 0.90 edges: %.3f + %.3f - %.3f\n"%(Mf_median, Mf_90max - Mf_median, Mf_median - Mf_90min))

  # Now do chif

  # Outputting the median and maP values and 68%, 90%, and 95% confidence intervals for the original posterior
                                                                
  chif_median = np.median(chif)
  chif_maP = aux.cen_from_edg(chif_bins)[np.nanargmax(P_chif)]

  chif_conf = aux.confidence_1d(P_chif, chif_bins)

  #chif_50min, chif_50max = chif_conf.edges(0.5)                                                                                                                    
                                                                                                                                                  
  chif_68min, chif_68max = chif_conf.edges(0.68)
  chif_90min, chif_90max = chif_conf.edges(0.90)
  chif_95min, chif_95max = chif_conf.edges(0.95)

  print("Original values:")
  print("chif median: %.3f"%(chif_median))
  print("chif maP: %.3f"%(chif_maP))
  #print("chif 0.50 region: %.3f, %.3f"%(chif_50min, chif_50max))                                                                                     
                                                                                                                           
  print("chif 0.68 region: %.3f, %.3f"%(chif_68min, chif_68max))
  print("chif 0.90 region: %.3f, %.3f"%(chif_90min, chif_90max))
  print("chif 0.95 region: %.3f, %.3f\n"%(chif_95min, chif_95max))

  print("chif median \pm 0.90 edges: %.3f + %.3f - %.3f\n"%(chif_median, chif_90max - chif_median, chif_median - chif_90min))

  return Mf, chif
