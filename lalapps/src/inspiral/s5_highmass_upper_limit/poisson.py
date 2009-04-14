#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
from pylab import *
from scipy.special import gamma,gammainc

xmax = 10.**3.
xmin = 10.**-2.
x = (xmax/xmin)**arange(0.,1.,0.001)*xmin

ymax = 10.**2.
ymin = 10**-2.
y = (ymax/ymin)**arange(0.,1.,0.001)*ymin
X,Y = meshgrid(x,y)
Z = X*0.

s1u = {}
s1d = {}
s2u = {}
s2d = {}
s3u = {}
s3d = {}
for s,p in zip([s1u, s2u, s3u],[0.158, 0.022, 0.001]):
  s["x"] = []
  s["y"] = []
  s["p"] = p

for s,p in zip([s1d, s2d, s3d],[0.158, 0.022, 0.001]):
  s["x"] = []
  s["y"] = []
  s["p"] = p

for idx in range(len(X[0])):
  u = X[:,idx]
  v = Y[:,idx]
  z = (v * (1./u)**(v) * exp(-1. / u) / gamma(v + 1.))
  z0 = exp(-1. / u[0])


  czd = gammainc(v + 1., 1. / u[0])
  czu = 1. - czd

  czd = [max(thisc, 10**-6) for thisc in czd]
  czu = [max(thisc, 10**-6) for thisc in czu]

  for s in [s1d, s2d, s3d]:
    if czd[0] >= 1. - s["p"]:
      s["x"].append(u[0])
      s["y"].append([thisv for thisz,thisv in zip(czd,v) if thisz >= 1. - s["p"]][-1])

  for s in [s1u, s2u, s3u]:
    if czu[-1] >= 1. - s["p"]:
      s["x"].append(u[0])
      s["y"].append([thisv for thisz,thisv in zip(czu,v) if thisz >= 1. - s["p"]][0])

  z /= max(z)
  Z[:,idx] = z


colors = []
dc = 1./200.
for idx in arange(0.,.5+dc,dc):
  colors.append((1.-idx,1.-idx,1.-idx))

figure()
hold(True)
subplot(111)
loglog(1./y,y,'k--')
for s in [s1u, s1d, s2u, s2d, s3u, s3d]:
  loglog(s["x"],s["y"],'k')
contourf(x,y,Z,len(colors),colors=colors)
ylim((0.9,ymax))
ylabel('Cumulative Number',size='x-large')
xlabel('IFAN',size='x-large')
savefig('Poisson_pdf_loglog.png')

X = load("val_zerolag.txt")
FAN = X[:,1]/25.
num = arange(len(FAN))+1
num = list(num)
num.reverse()

figure()
hold(True)
subplot(111)
loglog(1./y,y,'k--')
for s in [s1u, s1d, s2u, s2d, s3u, s3d]:
  loglog(s["x"],s["y"],'k')
loglog(1./FAN, num, '2', markersize=10, markeredgewidth=2)
contourf(x,y,Z,len(colors),colors=colors)
ylim((0.9,ymax))
ylabel('Cumulative Number',size='x-large')
xlabel(r'IFAN(MVSC)',size='x-large')
subplots_adjust(bottom=.125)
savefig('Poisson_pdf_loglog_MVSC.png')


X = load("val_zerolag_effsnr.txt")
FAN = X[:,1]/25.
num = arange(len(FAN))+1
num = list(num)
num.reverse()

figure()
hold(True)
subplot(111)
loglog(1./y,y,'k--')
for s in [s1u, s1d, s2u, s2d, s3u, s3d]:
  loglog(s["x"],s["y"],'k')
loglog(1./FAN, num, '2', markersize=10, markeredgewidth=2)
contourf(x,y,Z,len(colors),colors=colors)
ylim((0.9,ymax))
ylabel('Cumulative Number',size='x-large')
xlabel(r'IFAN($\rho_{\rm eff}$)',size='x-large')
subplots_adjust(bottom=.125)
savefig('Poisson_pdf_loglog_effsnr.png')


