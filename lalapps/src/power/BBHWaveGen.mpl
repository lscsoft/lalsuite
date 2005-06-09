interface(echo=0,quiet=true,prompt="->"):
##############################################################################
#
#  Warren's Binary Blackhole Waveform Generator
#
#  Author: Warren G. Anderson  (warren@phys.utb.edu)
#
#  Last Modified : Mon Jan 12 02:00 PM 2004 M
#
##############################################################################
# The idea of this program is to allow one to generate BBH kludge waveforms in
# a semi-automated way. The calculation of the waveform follows the paper of
# Flanagan and Hughes (F&H) to get h_+ and h_x for a binary black hole system
# that includes a 2PN inspiral, an ad-hoc merger which has characteristics
# similar to simulations that currently exist, and a ringdown. The masses and
# final angular momentum of the binary system will be parameters we select. At
# the end of the day, I hope to code some C functions that I can use to
# generate these signals on the fly.

# Constants and Units
kpc:=3.086e19*m:
c:=3.00e8*m/s:
T_sol:=4.925e-6*s:

# Read parameter file. If read fails, print error message and quit.

if assigned(filename)=false then
   filename:="BBHWaveGen.in"
fi:
try
   read filename:
catch:
   printf("ERROR: file %s not readable. ",filename):
   printf("Please check for existence and permissions.\n"):
   quit:
   error:
end try:
printf("Proceeding with input file %s.\n",filename):
   
# If a parameter do not exist or look wrong, alert user and quit.

if (not assigned(f_low) ) then
   printf("No value found for parameter f_low in file %s.\n",filename):
   quit:
fi:

if (not (type(f_low*s,numeric) and (f_low*s > 0)) ) then
   printf("f_low must be a positive number divided by the symbol s - "):
   printf("got %a instead.\n",f_low):
   quit:
fi:

if (not assigned(m1) ) then
   printf("No value found for parameter m1 in file %s.\n",filename):
   quit:
fi:

if (not (type(m1/M_sol,numeric) and (m1/M_sol > 1)) ) then
   printf("m1 must be a positive number multiplied by the symbol M_sol - "):
   printf("got %a instead.\n",m1):
   quit:
fi:

if (not assigned(m2) ) then
   printf("No value found for parameter m2 in file %s.\n",filename):
   quit:
fi:

if (not (type(m2/M_sol,numeric) and (m2/M_sol > 1)) ) then
   printf("m2 must be a positive number multiplied by the symbol M_sol - "):
   printf("got %a instead.\n",m2):
   quit:
fi:

if (not assigned(a) ) then
   printf("No value found for parameter a in file %s.\n",filename):
   quit:
fi:

if (not (type(a,numeric) and ((a>=0) and (a<=1))) ) then
   printf("a must be a number between 0 and 1 - "):
   printf("got %a instead.\n",a):
   quit:
fi:

if (not assigned(R) ) then
   printf("No value found for parameter R in file %s.\n",filename):
   quit:
fi:

if (not (type(R/m,numeric) and (R/m>0)) ) then
   printf("R must be a number multiplied by the symbol kpc or the symbol m - "):
   printf("got %a instead.\n",R):
   quit:
fi:

if (not assigned(eps) ) then
   printf("No value found for parameter eps in file %s.\n",filename):
   quit:
fi:

if(not (type(eps,numeric) and ((eps>0) and (eps<1))) ) then
   printf("eps must be a number between 0 and 1 - "):
   printf("got %a instead.\n",eps):
   quit:
fi:

if (not assigned(E_ratio) ) then
   printf("No value found for parameter E_ratio in file %s.\n",filename):
   quit:
fi:

if(not (type(E_ratio,numeric) and (E_ratio>0)) ) then
   printf("E_ratio must be a positive number - "):
   printf("got %a instead.\n",E_ratio):
   quit:
fi:

# Calculate total system mass and eta (reduced mass / total mass)
M_tot:=m1+m2:
eta:=m1*m2/M_tot^2:

# Calculate Inspiral Waveform

# First phase and frequency. We use the 2PN inspiral waveform given in the
# GRASP manual in section 6.4. Note that I am using gravitational wave
# frequency and phase (rather than the orbital values given in GRASP).
Theta:=eta*M_sol*(tc-t)/(5*T_sol*M_tot):
f_chirp:=expand((Theta^(-3/8)+((743/2688)+(11/32)*eta)*Theta^(-5/8)
   -(3*Pi/10)*Theta^(-3/4)+((1855099/14450688)+(56975/258048)*eta
   +(371/2048)*eta^2)*Theta^(-7/8))*M_sol/(8*Pi*T_sol*M_tot)):

# To set tc I use the estimate of the frequency at ISCO given by
# Flanagan and Hughes (eq (3.3)).
f_ISCO:=205/s*20*M_sol/M_tot:

# Now calculate the coalescence time so that the ISCO will be reached at
# t=0.
tc_list:=[solve(subs({t=0},f_chirp=f_ISCO),tc)]:

# The solver finds several solutions, so we need to check which is the
# actual value we want. We require it to be a value that has only
# positive frequencies before ISCO and for which the derivative of the
# frequency at the ISCO with respect to time is positive (frequency
# still increasing at the ISCO). I have nerver seen more than one of the
# three solutions satisfy these criteria, but I have not proved that
# this fortunate circumstance is necessary.
for i from 1 to nops(tc_list) do
   df_dtc:=diff(f_chirp,t):
   df_dtc:=evalf(subs(s=1,subs({tc=tc_list[i]},t=0,df_dtc))):
   f_eps:=evalf(subs(s=1,subs(tc=tc_list[i],t=-1e-4,f_chirp))):
   if ((df_dtc > 0) and f_eps > 0) then
      tc:=tc_list[i]:
      break:
   fi
od:i:='i':

# Now I can calculate when the signal enters the initial LIGO band, at a
# user specified frequency.
t0_list:=[solve(expand(f_chirp=40/s),t)]:
# The solver finds several solutions, so we need to check which is the
# actual one we want. Since we are setting t=0 at the ISCO, the start
# time must be negative. Fortunately, again, I've never seen more than
# one negative solution, so I can use that one.
for i from 1 to nops(t0_list) do
   if subs(s=1,t0_list[i])<0 then
      t0:=t0_list[i]:
      break:
   fi
od: i:='i':
# Integrate the frequency to find the phase evolution of the chirp.
phi_chirp:=2*Pi*expand(int(f_chirp,t)):
# Find the phase offset that would make the phase have value 0 at t=0
# and add it to the phase.
phi_chirp_0:=evalf(subs(t=0,phi_chirp)/Pi)*Pi:
phi_chirp:=phi_chirp-phi_chirp_0:
#
# h_+ and h_x for the Chirp
# Calculate h_plus and h_cross as per equations (6.6.3) and (6.6.4) of
# the GRASP manual. Note that I have assumed an optimal inclination
# angle. Also, I must divide the frequency by 2 in order to get the
# orbital frequency, which is what is used in these formulae.
A_chirp:=(T_sol*c/R)*(M_tot/M_sol)*
   expand(4*eta*(2*Pi*T_sol*M_tot*(f_chirp/2)/M_sol)^(2/3)):
hc_plus:=A_chirp*cos(phi_chirp):
hc_cross:=A_chirp*sin(phi_chirp):
# Ringdown
# Now for the ringdown portion.
# Parameter Values
# The following parameter values are extracted from the LAL Software
# Document ring package.
Q:=2*(1-a)^(-9/10):
f_qnr:=3.2e4*(1-0.63*(1-a)^(3/10))*(M_sol/M_tot)/s:
A_qnr:=2.415e-21*(1e3*kpc/R)*(M_tot/M_sol)*sqrt(eps/0.01)
   /(sqrt(2*Pi*Q)*sqrt(1-0.63*(1-a)^(3/10)))*exp(-Pi*f_qnr*t/Q):
#
# h_plus and h_cross
# These are taken from the LSD.
hr_plus:=A_qnr*cos(2*Pi*f_qnr*t+phi_qnr_0):
hr_cross:=A_qnr*sin(2*Pi*f_qnr*t+phi_qnr_0):
# How long does this last? We allow five e-foldings.
dt_qnr:=5*Q/(Pi*f_qnr):
# Energy radiated by Ringdown
# This will be useful to know shortly. It is simply proportional the
# integral of the square of dh/dt.
Er_plus:=int(diff(subs(sec=1,hr_plus),t)^2,t=0..dt_qnr):
Er_cross:=int(diff(subs(sec=1,hr_cross),t)^2,t=0..dt_qnr):
#
# Merger
# This is my cooked up merger waveform. I am trying to make it
# consistent with the predictions of F&H. I therefore require it to have
# ~3 times the energy of the ringdown, to last for 50M seconds, and to
# have frequency components between f_ISCO and f_qnr. I will take the
# general form of the waveform to be
# h(t):=h_avg_merge(t)*exp(2*Pi*I*f(t)).
# Merger Duration
# From just before (3.26) of F&H.
dt_merge:=50*(M_tot/M_sol)*T_sol:
#
# Frequency Evolution
# I don't know what the frequency evolution should be, but I know it
# should be in the band from
# f_ISCO < f < f_qnr and I would like the frequency and it's derivative
# with respect to time to be continuous at both ends. This is actually a
# set of four conditions, so I need four parameters in my curve with
# which to fit them, so I fit a cubic, which is the easiest 4 parameter
# curve I can think of.
f_merge:=f0+f1*t+f2*t^2+f3*t^3:
# Fit the frequency and its derivative from the left.
f0:=subs(t=0,f_chirp):
f1:=solve(subs(t=0,diff(f_merge,t))=subs(t=0,diff(f_chirp,t)),f1):
# Fit the frequency and its derivative from the right.
eq1:=subs(t=dt_merge,f_merge)=f_qnr:
eq2:=subs(t=dt_merge,diff(f_merge,t))=0:
f_co:=convert(solve({eq1,eq2},{f3,f2}),list):
assign(eval(f_co[1])):
assign(eval(f_co[2])):
#
# Phase Evolution and final phase
# I now integrate the the frequency to get the phase evolution and the
# final phase.
phi_merge:=2*Pi*int(f_merge,t):
phi_qnr_0:=evalf(subs(t=dt_merge,phi_merge)/Pi)*Pi:
#
# Average Amplitude Evolution
# Again, lacking any better knowledge of how to do this, I assume the
# angle averaged rms amplitude to match h_avg_chirp and its time
# deriviative at t=0 and h_qnr(phi_merge_final) and its time derivative
# at t=T_merge. I also want it to have 3 times the energy of the
# ringdown in the merger waveform, which puts a further restriction on
# the amplitude evolution. Thus, there are five constraints, and I need
# a five parameter function. The simplest five parameter function I can
# think of is a quartic, so that's what I'll use.
A_merge:=A0+A1*t+A2*t^2+A3*t^3+A4*t^4:
# Match from the left.
A0:=evalf(subs(t=0,A_chirp)):
A1:=evalf(subs(t=0,diff(A_chirp,t))):
# Match from the right.
A2:=solve(subs(t=dt_merge,A_merge)=subs(t=0,A_qnr),A2):
A3:=solve(subs(t=dt_merge,diff(A_merge,t))=subs(t=0,diff(A_qnr,t)),A3):
#
# Waveform (almost)
# Calculate the waveforms apart from the unkown constant A4.
hm_plus:=A_merge*cos(phi_merge):
hm_cross:=A_merge*sin(phi_merge):
#
#
# Energy Condition
# I still have one coefficient of h_avg_merge to fix. I do this by
# requiring the energy to be a user-defined multiple of the qnr energy.
# Flanagan and Hughes have a factor of about 3 for an a=0.98 final BH.
Em_plus:=evalf(Int(subs({s=1,h4=3.9},diff(hm_plus,t)^2),t=0..dt_merge/s)):
Em_cross:=evalf(Int(subs({s=1,h4=3.9},diff(hm_cross,t)^2),t=0..dt_merge/s)):
# These need to be set to three times the energies in the qnr sections.
# We therefore solve:
A4:=fsolve(subs(s=1,Em_plus=E_ratio*Er_plus),A4):
# Now check that this is ok for the other polarization.
evalf(subs(s=1,Em_cross/Er_cross)):
#
# Merger Amplitude Evolution and Waveform
#
# Combined Waveform for BH-BH system
# Now to put it all together.The chirp starts in the infinite past
# (although we cut it off at a user defined lower frequency) and ends at
# t=0 and with phi=0, as I have already arranged. The merger starts at
# phi=0 and at t=0 and goes to t=dt_merge, where it ends with
# phi=phi_qnr_0. The ringdown goes from t=dt_merge and phi=phi_qnr_0 to
# the infinite future, although again we have imposed a cutoff at five
# e-foldings of the ringdown..
t_end:=evalf(dt_merge+dt_qnr):
hc_plus:=hc_plus*Heaviside(-t):
hc_cross:=hc_cross*Heaviside(-t):
hm_plus:=hm_plus*Heaviside(t)*Heaviside(dt_merge/s-t):
hm_cross:=hm_cross*Heaviside(t)*Heaviside(dt_merge/s-t):
hr_plus:=subs(t=t-dt_merge,hr_plus)*Heaviside(t-dt_merge/s):
hr_cross:=subs(t=t-dt_merge,hr_cross)*Heaviside(t-dt_merge/s):
h_plus:=hc_plus+hm_plus+hr_plus:
h_cross:=hc_cross+hm_cross+hr_cross:
#
# Output to File
# I need to write out the data to make it easier to use. I choose my
# sampling frequency so that the Nyquist frequency is comfortably higher
# than the band where my signal has significant power. I choose the
# sampling frequency to be a power of 2 because my conciousness has been
# permanently warped by doing FFT's and it seems most natural to do so
# :).
ofilename:=cat(filename,".dat"):
try of:=fopen(ofilename,WRITE)
catch:
   printf("ERROR: file %s cannot be opened. ",ofilename):
   printf("Please check permissions.\n"):
   quit:
   error:
end try:
# I can now calculate the number of time samples I will use.
f_samp:=2^(ceil(log[2](f_qnr*s))+1)/s:
N:=round(f_samp*(t_end-t0)):
s:=1:
for i from 0 to N-1 do
   ti:=evalf(i/f_samp):
   h_plusi:=evalc(Re(evalf(subs(t=t0+ti,h_plus)))):
   h_crossi:=evalc(Re(evalf(subs(t=t0+ti,h_cross)))):
   fprintf(of,"%e %e %e\n",ti,h_plusi,h_crossi):
od: i:='i':
fclose(of):
printf("Output printed to file %s.\n",ofilename):
quit:
