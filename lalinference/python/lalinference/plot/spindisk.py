# Ben Farr 2016

import numpy as np
import os
from matplotlib.lines import Line2D

__all__=('make_disk_plot',)

def make_disk_plot(post,outpath=None):
  from matplotlib import pyplot as plt
  from matplotlib import rc
  rc('text', usetex=False)
  rc('font', family='lmodern')

  small, big = 10, 12
  rc('axes', labelsize=big)
  rc('font', size=big)
  rc('legend', fontsize=small)
  rc('xtick', labelsize=small)
  rc('ytick', labelsize=small)
  rc('text.latex', unicode=True)


  try:
    import corner
  except ImportError:
    print("cannot import corner. Won't plot spin disk")
    return None

  a1='a1'
  tilt1='tilt1'
  a2='a2'
  tilt2='tilt2'

  names=post.names
  
  if not set([a1,a2,tilt1,tilt2]).issubset(names):
    print("Cannot plot spin disk plot. Not all required spin parameters exist in the posterior file. Skipping...\n")
    return None

  # Spin disk plot
  fig, axs = plt.subplots(1, 2, sharey=True, figsize=(4, 4))

  Na, Nt = 20, 30
  xticks = [0., .25, .5, .75, 1.]

  vmin, vmax = 0., 0.
  for a, tilt in zip([a1, a2], [tilt1, tilt2]):
      asamps=(post[a].samples).flatten()
      tsamps=(post[tilt].samples).flatten()
      
      H, _, _ = np.histogram2d(asamps, np.cos(tsamps), range=[[0, 1], [-1, 1]], bins=(Na, Nt), normed=False)
      H /= len(asamps)
      vmax = H.max() if H.max() > vmax else vmax

  for ax, a, tilt, flip in zip(axs, [a1, a2], [tilt1, tilt2], [True, False]):
      plt.sca(ax)
      H, rs, costs = np.histogram2d(asamps, np.cos(tsamps), range=[[0, 1], [-1, 1]], bins=(Na, Nt), normed=False)
      H /= len(asamps)
      COSTS, RS = np.meshgrid(costs, rs)
      X = RS * np.sin(np.arccos(COSTS))
      Y = RS * COSTS

      HS = np.column_stack((X.flatten(), Y.flatten()))

      XS = np.reshape(HS[:,0], (Na+1,Nt+1))
      YS = np.reshape(HS[:,1], (Na+1,Nt+1))

      plt.pcolormesh(XS, YS, H, vmin=vmin, vmax=vmax, edgecolor='face', cmap='Greys')

      ax.set_ylim((-1., 1.))
      ax.set_xlim((0., 1.))
      if flip:
          ax.set_xticks(xticks[1:])
          ax.invert_xaxis()
      else:
          ax.set_xticks(xticks)
          ax.yaxis.tick_right()
          ax.yaxis.set_label_position("right")

  axs[0].set_xlabel(r'$|\mathbf{S_1} \times \mathbf{\hat{L}}|$')
  axs[1].set_xlabel(r'$|\mathbf{S_2} \times \mathbf{\hat{L}}|$')
  axs[0].set_ylabel('$\mathbf{S_1}\cdot\mathbf{\hat{L}}$')
  axs[1].set_ylabel('$\mathbf{S_2}\cdot\mathbf{\hat{L}}$')

  fig.subplots_adjust(wspace=0.04)
  cax = fig.add_axes([0.06, -0.075, 0.9, 0.05])

  cbar = plt.colorbar(orientation='horizontal', cax=cax)
  cbar.formatter.set_powerlimits((-1, 1))
  cbar.update_ticks()
  cbar.set_label('posterior probability')
  cbar.solids.set_edgecolor("face")
  plt.savefig(os.path.join(outpath,"comp_spin_pos.png"), bbox_inches='tight')
  plt.clf()

