/************************** <lalVerbatim file="GeneratePPNInspiralCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GeneratePPNInspiral.c}}
\label{ss:GeneratePPNInspiral.c}

Computes a parametrized post-Newtonian inspiral waveform.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GeneratePPNInspiralCP}
\index{\texttt{LALGeneratePPNInspiral()}}

\subsubsection*{Description}

This function computes an inspiral waveform using the parameters in
\verb@*params@, storing the result in \verb@*output@.

In the \verb@*params@ structure, the routine uses all the ``input''
fields specified in \verb@GeneratePPNInspiral.h@, and sets all of the
``output'' fields.  If \verb@params->ppn=NULL@, a normal
post${}^2$-Newtonian waveform is generated; i.e.\ $p_0=1$, $p_1=0$,
$p_2=1$, $p_3=1$, $p_4=1$, $p_{5+}=0$.

In the \verb@*output@ structure, the fields \verb@output->a@ and
\verb@output->phi@ must exist, since the routine uses
\verb@output->a->deltaT@ and \verb@output->phi->deltaT@ to determine
the sampling rate.  Their values \emph{must} be the same.  The fields
\verb@output->a->epoch@ and \verb@output->phi->epoch@ are arbitrary,
but must also be the same.  If \verb@output->a->data@ and
\verb@output->phi->data@ are both set, \emph{and}
\verb@output->a->data->vectorLength@=2, \emph{and}
\verb@output->phi->data->vectorLength@=1, \emph{and}
\verb@output->a->data->length@=\verb@output->phi->data->length@, then
these arrays will be filled with the waveform (and padded with zeros
if the waveform terminates before the end of the array).  Otherwise,
\verb@output->a->data@ and \verb@output->phi->data@ will be destroyed
(if necessary) and reallocated to a length appropriate to the waveform
generated.  No other fields in \verb@*output@ are examined or changed.

\subsubsection*{Algorithm}

This function is a fairly straightforward calculation of
Eqs.~\ref{eq:ppn-freq}--\ref{eq:ppn-across} in
\verb@GeneratePPNInspiral.h@.  However, there are some nontrivial
issues involved, which are discussed in some depth in Secs.~6.4, 6.6,
and~6.9.2 of~\cite{GRASP_1.9.8:2000}.  What follows is a brief
discussion of these issues and how this routine deals with them.

\vspace{1ex}\noindent\textbf{Computing the start time:}

When building a waveform for data analysis, one would generally like
to start the waveform at some well-defined frequency where it first
enters the band of interest; one then defines the start time of the
integration by inverting Eq.~\ref{eq:ppn-freq} if
\verb@GeneratePPNInspiral.h@.  The current algorithm follows this
standard approach by requiring the calling routine to specify
\verb@params->fStartIn@, which is then inverted to find
$\Theta_\mathrm{start}$.  This inversion is in fact the most
algorithmically complicated part of the routine, so we will discuss it
in depth.

To help clarify the problem, let us rewrite the equation in
dimensionless parameters $y=8\pi fT_\odot m_\mathrm{tot}/M_\odot$ and
$x=\Theta^{-1/8}$:
\begin{equation}
\label{eq:ppn-fnorm}
y = \sum_{k=0}^{5} C_k x^{k+3} \; ,
\end{equation}
where:
\begin{eqnarray}
C_0 & = & p_0 \;,\nonumber\\
C_1 & = & p_1 \;,\nonumber\\
C_2 & = & p_2\left(\frac{743}{2688}+\frac{11}{32}\eta\right) \;,\nonumber\\
C_3 & = & p_3\frac{3\pi}{10} \;,\nonumber\\
C_4 & = & p_4\left(\frac{1855099}{14450688}+\frac{56975}{258048}\eta+
		\frac{371}{2048}\eta^2\right) \;,\nonumber\\
C_5 & = & p_5\left(\frac{7729}{21504}+\frac{3}{256}\eta\right)\pi \;.\nonumber
\end{eqnarray}
We note that $x$ is a time parameter mapping the range
$t=-\infty\rightarrow t_c$ to $x=0\rightarrow\infty$.

In a normal post-Newtonian expansion it is possible to characterize
the general behaviour of this equation quite accurately, since the
values of $p_k$ are known and since $\eta$ varies only over the range
$[0,1/4]$.  In a parametrized post-Newtonian expansion, however, even
the relative orders of magnitude of the coefficients can vary
significantly, making a robust generic root finder impractical.
However, we are saved by the fact that we can restrict our search to
domains where the post-Newtonian expansion is a valid approximation.
We define the post-Newtonian expansion \emph{not} to be valid if
\emph{any} of the following conditions occur:
\begin{enumerate}
\item A higher-order term in the frequency expansion becomes larger in
magnitude than the leading (lowest-order nonzero) term.
\item The inferred orbital radius, approximated by
$r\sim4m_\mathrm{tot}\Theta^{1/4}$, drops below $2m_\mathrm{tot}$;
i.e.\ $\Theta<1/16$ or $x>1/\sqrt{2}$.
\item The frequency evolution becomes non-monotonic.
\end{enumerate}
We can further require as a matter of convention that the lowest-order
nonzero coefficient in the frequency expansion be positive; this is
simply a sign convention stating that the frequency of a system be
positive at large radii.

The first two conditions above allow us to set firm limits on the
range of the initial $x_\mathrm{start}$.  Let $C_j$ be the
lowest-order nonzero coefficient; then for every nonzero $C_{k>j}$
we can define a point $x_k=|C_j/C_k|^{1/(k-j)}$ where that term
exceeds the leading-order term in magnitude.  We can therefore limit
the range of $x$ to values less than $x_\mathrm{max}$, which is the
minimum of $1/\sqrt{2}$ and all $x_k$.  We note that even if we were
to extend the post-Newtonian expansion in Eq.~\ref{eq:ppn-fnorm} to an
infinite number of terms, this definition of $x_\mathrm{max}$ implies
that the frequency is guaranteed to be monotonic up to
$x_\mathrm{max}(5-\sqrt{7})/6$, and positive up to $x_\mathrm{max}/2$.
Thus we can confidently begin our search for $x_\mathrm{start}$ in the
domain $(0,0.39x_\mathrm{max})$, where the leading-order term
dominates, and end it if we ever exceed $x_\mathrm{max}$.

We therefore bracket our value of $x_\mathrm{start}$ as follows: We
start with an initial guess
$x_\mathrm{guess}=(y_\mathrm{start}/C_j)^{1/(j+3)}$, or
$0.39x_\mathrm{max}$, whichever is less.  If
$y(x_\mathrm{guess})<y_\mathrm{start}$, we iteratively decrease $x$ by
factors of 0.95 until $y(x)>y_\mathrm{start}$; this is guaranteed to
occur within a few iterations, since we are moving into a regime where
the leading-order behaviour dominates more and more.  If
$y(x_\mathrm{guess})>y_\mathrm{start}$, we iteratively increase $x$ by
factors of 1.05 until $y(x)>y_\mathrm{start}$, or until
$x>x_\mathrm{max}$; this is also guaranteed to occur quickly because,
in the worst case, it only takes about 20 iterations to step from
$0.39x_\mathrm{max}$ to $x_\mathrm{max}$, and if $x_\mathrm{guess}$
were much lower than $0.39x_\mathrm{max}$ it would have been a pretty
good guess to begin with.  If at any point while increasing $x$ we
find that $y$ is decreasing, we determine that the starting frequency
is already in a regime where the post-Newtonian approximation is
invalid, and we return an error.  Otherwise, once we have bracketed
the value of $x_\mathrm{start}$, we use \verb@LALSBisectionFindRoot()@
to pin down the value to an accuracy of a part in $10^6$.

\vspace{1ex}\noindent\textbf{Computing the phase and amplitudes:}

Once we have $x_\mathrm{start}$, we can find $\Theta_\mathrm{start}$,
and begin incrementing it; at each timestep we compute $x$ and hence
$f$, $\phi$, $A_+$, and $A_\times$ according to
Eqs.~\ref{eq:ppn-freq}, \ref{eq:ppn-phi}, \ref{eq:ppn-aplus},
and~\ref{eq:ppn-across}.  If valid output arrays were provided in
\verb@*output@, these will be filled; otherwise the routine
progressively creates a list of length-1024 arrays and fills them.
The process stops when any of the following occurs:
\begin{enumerate}
\item The frequency exceeds the requested termination frequency.
\item The number of steps reaches the length of the provided arrays
(if valid arrays were provided), or the suggested maximum length in
\verb@*params@, whichever is less.
\item The frequency is no longer increasing.
\item The parameter $x>x_\mathrm{max}$.
\item We run out of memory.
\end{enumerate}
In the last case an error is returned; otherwise the waveform is
deemed ``complete''.  If arrays were provided any remaining space is
padded with zeros; otherwise, output arrays are created of the
appropriate length and are filled with the data.

Internally, the routine keeps a list of all coefficients, as well as a
list of booleans indicating which terms are nonzero.  The latter
allows the code to avoid lengthy floating-point operations (especially
the logarithm in the post${}^{5/2}$-Newtonian phase term) when these
are not required.

\vspace{1ex}\noindent\textbf{Warnings and suggestions:}

If no post-Newtonian parameters are provided (i.e.\
\verb@params->ppn=NULL@), we generate a post${}^2$-Newtonian waveform,
\emph{not} a post${}^{5/2}$-Newtonian waveform.  This is done not only
for computationally efficiency, but also because the accuracy and
reliability of the post${}^{5/2}$-Newtonian waveform is actually
worse.  You can of course specify a post${}^{5/2}$-Newtonian waveform
with an appropriate assignment of \verb@params->ppn@, but you do so at
your own risk!

This routine also performs no sanity checking on the requested
sampling interval $\Delta
t=$\verb@output->a->deltaT@=\verb@output->phi->deltaT@, because this
depends very much on how one intends to use the generated waveform.
If you plan to generate actual wave functions $h_{+,\times}(t)$ at the
same sample rate, then you will generally want a sampling interval
$\Delta t<1/2f_\mathrm{max}$; you can enforce this by specifying a
suitable \verb@params->fStopIn@.

However, many routines (such as those in \verb@SimulateCoherentGW.h@)
generate actual wave functions by linear interpolation of the
amplitude and phase data, which then need only be sampled on
timescales $\sim\dot{f}^{-1/2}$ rather than $\sim f^{-1}$.  More
precisely, we would like our interpolated phase to differ from the
actual phase by no more than some specified amount, say $\pi$ radians.
The largest deviation from linear phase evolution will typically be on
the order of $\Delta\phi\approx(1/2)\ddot{\phi}(\Delta
t/2)^2\approx(\pi/4)\Delta f\Delta t$, where $\Delta f$ is the
frequency shift over the timestep.  Thus in general we would like to
have
$$
\Delta f \Delta t \lessim 4
$$
for our linear interpolation to be valid.  This routine helps out by
setting the output parameter field \verb@params->dfdt@ equal to the
maximum value of $\Delta f\Delta t$ encountered during the
integration.

\subsubsection*{Uses}
\begin{verbatim}
LALSCreateVectorSequence()
LALSDestroyVectorSequence()
LALSBisectionFindRoot()
LALWarning()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GeneratePPNInspiralCV}}

******************************************************* </lalLaTeX> */

/*********************************************************************
 * PREAMBLE                                                          *
 *********************************************************************/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID( GENERATEPPNINSPIRALC, "$Id$" );

/* Define some constants used in this module. */
#define MAXORDER 6        /* Maximum number of N and PN terms */
#define BUFFSIZE 1024     /* Number of timesteps buffered */
#define ACCURACY (1.0e-6) /* Accuracy of root finder */
#define TWOTHIRDS (0.6666666667) /* 2/3 */

/* A macro to computing the (normalized) frequency.  It appears in
   many places, including in the main loop, and I don't want the
   overhead of a function call.  The following variables are required
   to be defined and set outside of the macro:

   REAL4 c0, c1, c2, c3, c4, c5;   PN frequency coefficients
   BOOLEAN b0, b1, b2, b3, b4, b5; whether to include each PN term

   The following variables must be defined outside the macro, but are
   set inside it:

   REAL4 x2, x3;  the square and cube of the input x */
#define FREQ( f, x )                                                 \
do {                                                                 \
  x2 = (x)*(x);                                                      \
  x3 = x2*(x);                                                       \
  (f) = 0;                                                           \
  if ( b0 )                                                          \
    (f) += c0;                                                       \
  if ( b1 )                                                          \
    (f) += c1*(x);                                                   \
  if ( b2 )                                                          \
    (f) += c2*x2;                                                    \
  if ( b3 )                                                          \
    (f) += c3*x3;                                                    \
  if ( b4 )                                                          \
    (f) += c4*x3*(x);                                                \
  if ( b5 )                                                          \
    (f) += c5*x3*x2;                                                 \
  (f) *= x3;                                                         \
} while (0)

/* Definition of a data structure used by FreqDiff() below. */
typedef struct tagFreqDiffParamStruc {
  REAL4 *c;   /* PN coefficients of frequency series */
  BOOLEAN *b; /* whether to include each PN term */
  REAL4 y0;   /* normalized frequency being sought */
} FreqDiffParamStruc;

/* A function to compute the difference between the current and
   requested normalized frequency, used by the root bisector. */
static void
FreqDiff( LALStatus *stat, REAL4 *y, REAL4 x, void *p )
{
  INT4 i;     /* index over PN coefficients */
  REAL4 f;    /* normalized frequency */
  REAL4 *c;   /* PN coefficients of frequency series */
  BOOLEAN *b; /* whether to include each PN term */

  INITSTATUS( stat, "Frequency", GENERATEPPNINSPIRALC );
  ASSERT( p, stat, 1, "Null pointer" );

  c = ( (FreqDiffParamStruc *)p )->c;
  b = ( (FreqDiffParamStruc *)p )->b;
  f = 0.0;
  for ( i = 0; i < MAXORDER; i++ )
    if ( b[i] )
      f += c[i]*pow( x, i + 3.0 );
  *y = f - ( (FreqDiffParamStruc *)p )->y0;
  RETURN( stat );
}

/* Definition of a data buffer list for storing the waveform. */
typedef struct tagPPNInspiralBuffer {
  REAL4 a[BUFFSIZE];                 /* amplitude data */
  REAL4 phi[2*BUFFSIZE];             /* phase data */
  struct tagPPNInspiralBuffer *next; /* next buffer in list */
} PPNInspiralBuffer;

/* Definition of a macro to free the tail of said list, from a given
   node onward. */
#define FREELIST( node )                                             \
do {                                                                 \
  PPNInspiralBuffer *herePtr = (node);                               \
  while ( herePtr ) {                                                \
    PPNInspiralBuffer *lastPtr = herePtr;                            \
    herePtr = herePtr->next;                                         \
    LALFree( lastPtr );                                              \
  }                                                                  \
} while (0)


/*********************************************************************
 * MAIN FUNCTION                                                     *
 *********************************************************************/

/* <lalVerbatim file="GeneratePPNInspiralCP"> */
void
LALGeneratePPNInspiral( LALStatus     *stat,
			CoherentGW    *output,
		        PPNParamStruc *params )
{ /* </lalVerbatim> */

  /* System-derived constants. */
  BOOLEAN b0, b1, b2, b3, b4, b5; /* whether each order is nonzero */
  BOOLEAN b[MAXORDER];            /* vector of above coefficients */
  REAL4 c0, c1, c2, c3, c4, c5;   /* PN frequency coefficients */
  REAL4 c[MAXORDER];              /* vector of above coefficients */
  REAL4 d0, d1, d2, d3, d4, d5;   /* PN phase coefficients */
  REAL4 p[MAXORDER];              /* PN parameter values */
  REAL4 mTot, mu;      /* total mass and reduced mass */
  REAL4 eta, etaInv;   /* mass ratio and its inverse */
  REAL4 phiC;          /* phase at coalescence */
  REAL4 cosI;          /* cosine of system inclination */
  REAL4 fFac;          /* SI normalization for f and t */
  REAL4 f2aFac;        /* factor multiplying f in amplitude function */
  REAL4 apFac, acFac;  /* extra factor in plus and cross amplitudes */

  /* Integration parameters. */
  BOOLEAN buffer;  /* whether we are buffering the data in a list */
  UINT4 i;         /* index over PN terms */
  UINT4 j;         /* index of leading nonzero PN term */
  UINT4 n, nMax;   /* index over timesteps, and its maximum + 1 */
  UINT4 nNext;     /* index where next buffer starts */
  REAL4 t, dt;     /* dimensionless time and increment */
  REAL4 x, xStart, xMax; /* x = t^(-1/8), and its maximum range */
  REAL4 y, yStart, yMax; /* normalized frequency and its range */
  REAL4 yOld, dyMax;     /* previous timestep y, and maximum y - yOld */
  REAL4 x2, x3;          /* x^2 and x^3 */
  REAL4 *a, *phi;        /* pointers to amplitude and phase data */
  PPNInspiralBuffer *head, *here; /* pointers to buffered data */

  INITSTATUS( stat, "LALGeneratePPNInspiral", GENERATEPPNINSPIRALC );
  ATTATCHSTATUSPTR( stat );

  /*******************************************************************
   * CHECK INPUT PARAMETERS                                          *
   *******************************************************************/

  /* Dumb initialization to shut gcc up. */
  head = here = NULL;
  b0 = b1 = b2 = b3 = b4 = b5 = 0;
  c0 = c1 = c2 = c3 = c4 = c5 = d0 = d1 = d2 = d3 = d4 = d5 = 0.0;

  /* Make sure parameter structures and their fields exist. */
  ASSERT( params, stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( output, stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );
  if ( !( ( output->a ) && ( output->phi ) ) ) {
    ABORT( stat, GENERATEPPNINSPIRALH_ESIG,
	   GENERATEPPNINSPIRALH_MSGESIG );
  }

  /* Make sure timing data in *output is consistent. */
  if ( ( output->a->deltaT != output->phi->deltaT ) ||
       ( output->a->epoch.gpsSeconds !=
	 output->phi->epoch.gpsSeconds ) ||
       ( output->a->epoch.gpsNanoSeconds !=
	 output->phi->epoch.gpsNanoSeconds ) ) {
    ABORT( stat, GENERATEPPNINSPIRALH_ETBAD,
	   GENERATEPPNINSPIRALH_MSGETBAD );
  }

  /* See whether output->a and output->phi are set up to store the
     data, or whether we'll have to buffer it. */
  buffer = !( output->a->data ) || !( output->phi->data ) ||
    !( output->a->data->data ) || !( output->phi->data->data ) ||
    ( output->a->data->length != output->phi->data->length ) ||
    ( output->a->data->vectorLength != 2 ) ||
    ( output->phi->data->vectorLength != 1 );

  /* Get PN parameters, if they are specified; otherwise use
     post2-Newtonian. */
  if ( params->ppn ) {
    ASSERT( params->ppn->data, stat, GENERATEPPNINSPIRALH_ENUL,
	    GENERATEPPNINSPIRALH_MSGENUL );
    j = params->ppn->length;
    if ( j > MAXORDER )
      j = MAXORDER;
    for ( i = 0; i < j; i++ )
      p[i] = params->ppn->data[i];
    for ( ; i < MAXORDER; i++ )
      p[i] = 0.0;
  } else {
    p[0] = 1.0;
    p[1] = 0.0;
    p[2] = 1.0;
    p[3] = 1.0;
    p[4] = 1.0;
    for ( i = 5; i < MAXORDER; i++ )
      p[i] = 0.0;
  }

  /*******************************************************************
   * COMPUTE SYSTEM PARAMETERS                                       *
   *******************************************************************/

  /* Compute parameters of the system. */
  mTot = params->mTot;
  ASSERT( mTot != 0.0, stat, GENERATEPPNINSPIRALH_EMBAD,
	  GENERATEPPNINSPIRALH_MSGEMBAD );
  eta = params->eta;
  ASSERT( eta != 0.0, stat, GENERATEPPNINSPIRALH_EMBAD,
	  GENERATEPPNINSPIRALH_MSGEMBAD );
  etaInv = 2.0 / eta;
  mu = eta*mTot;
  cosI = cos( params->inc );
  phiC = params->phi;

  /* Compute frequency, phase, and amplitude factors. */
  fFac = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
  dt = -output->a->deltaT * eta / ( 5.0*LAL_MTSUN_SI );
  ASSERT( dt > 0.0, stat, GENERATEPPNINSPIRALH_ETBAD,
	  GENERATEPPNINSPIRALH_MSGETBAD );
  f2aFac = LAL_TWOPI*LAL_MTSUN_SI*mTot*fFac;
  ASSERT( params->d != 0.0, stat, GENERATEPPNINSPIRALH_EDBAD,
	  GENERATEPPNINSPIRALH_MSGEDBAD );
  apFac = acFac = -2.0*mu*LAL_MRSUN_SI/params->d;
  apFac *= 1.0 + cosI*cosI;
  acFac *= 2.0*cosI;

  /* Compute PN expansion coefficients. */
  if ( p[0] != 0.0 ) {
    b0 = b[0] = 1;
    c0 = c[0] = p[0];
    d0 = p[0];
  } else
    b0 = b[0] = 0;
  if ( p[1] != 0.0 ) {
    b1 = b[1] = 1;
    c1 = c[1] = p[1];
    d1 = p[1]*5.0/4.0;
  } else
    b1 = b[1] = 0;
  if ( p[2] != 0.0 ) {
    b2 = b[2] = 1;
    c2 = c[2] = p[2]*( 743.0/2688.0 + eta*11.0/32.0 );
    d2 = p[2]*( 3715.0/8064.0 + eta*55.0/96.0 );
  } else
    b2 = b[2] = 0;
  if ( p[3] != 0.0 ) {
    b3 = b[3] = 1;
    c3 = c[3] = -p[3]*( 3.0*LAL_PI/10.0 );
    d3 = -p[3]*( 3.0*LAL_PI/4.0 );
  } else
    b3 = b[3] = 0;
  if ( p[4] != 0.0 ) {
    b4 = b[4] = 1;
    c4 = c[4] = p[4]*( 1855099.0/14450688.0 + eta*56975.0/258048.0 +
		       eta*eta*371.0/2048.0 );
    d4 = p[4]*( 9275495.0/14450688.0 + eta*284875.0/258048.0 +
		eta*eta*1855.0/2048.0 );
  } else
    b4 = b[4] = 0;
  if ( p[5] != 0.0 ) {
    b5 = b[5] = 1;
    c5 = c[5] = -p[5]*( 7729.0/21504.0 + eta*3.0/256.0 )*LAL_PI;
    d5 = -p[5]*( 38645.0/172032.0 + eta*15.0/2048.0 )*LAL_PI;
  } else
    b5 = b[5] = 0;

  /* Find the leading-order frequency term.  Note: This will work even
     if given negative (i.e. unphysical) values of eta. */
  for ( j = 0; ( j < MAXORDER ) && ( ( b[j] == 0 ) ||
				     ( c[j] == 0.0 ) ); j++ )
    ;
  if ( j == MAXORDER ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EPBAD,
	   GENERATEPPNINSPIRALH_MSGEPBAD );
  }

  /*******************************************************************
   * COMPUTE START TIME                                              *
   *******************************************************************/

  /* First, find the normalized start frequency, and the best guess as
     to the start time from the leading-order term.  We require the
     frequency to be increasing. */
  ASSERT( params->fStopIn > params->fStartIn, stat,
	  GENERATEPPNINSPIRALH_EFBAD, GENERATEPPNINSPIRALH_MSGEFBAD );
  yStart = params->fStartIn / fFac;
  yMax = params->fStopIn / fFac;
  if ( ( c[j]*fFac < 0.0 ) || ( yStart < 0.0 ) || (yMax < 0.0 ) ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EPBAD,
	   GENERATEPPNINSPIRALH_MSGEPBAD );
  }
  xStart = pow( yStart/c[j], 1.0/( j + 3.0 ) );
  xMax = LAL_SQRT1_2;

  /* The above is exact if the leading-order term is the only one in
     the expansion.  Check to see if there are any other terms. */
  for ( i = j + 1; ( i < MAXORDER ) && ( b[i] == 0 ); i++ )
    ;
  if ( i < MAXORDER ) {
    /* There are other terms, so we have to use bisection to find the
       start time. */
    REAL4 xLow, xHigh; /* ultimately these will bracket xStart */
    REAL4 yLow, yHigh; /* the normalized frequency at these times */

    /* If necessary, revise the estimate of the cutoff where we know
       the PN approximation goes bad, and revise our initial guess to
       lie well within the valid regime. */
    for ( i = j + 1; i < MAXORDER; i++ )
      if ( b[i] != 0 ) {
	x = pow( fabs( c[j]/c[i] ), 1.0/(REAL4)( i - j ) );
	if ( x < xMax )
	  xMax = x;
      }
    if ( xStart < 0.39*xMax )
      xStart = 0.39*xMax;

    /* If our frequency is too high, step backwards and/or forwards
       until we have bracketed the correct frequency. */
    xLow = xHigh = xStart;
    FREQ( yHigh, xStart );
    yLow = yHigh;
    while ( yLow > yStart ) {
      xHigh = xLow;
      yHigh = yLow;
      xLow *= 0.95;
      FREQ( yLow, xLow );
    }
    while ( yHigh < yStart ) {
      xLow = xHigh;
      yLow = yHigh;
      xHigh *= 1.05;
      FREQ( yHigh, xHigh );
      /* Check for PN breakdown. */
      if ( ( yHigh < yLow ) || ( xHigh > xMax ) ) {
	ABORT( stat, GENERATEPPNINSPIRALH_EFBAD,
	       GENERATEPPNINSPIRALH_MSGEFBAD );
      }
    }

    /* We may have gotten lucky and nailed the frequency right on.
       Otherwise, find xStart by root bisection. */
    if ( yLow == yStart )
      xStart = xLow;
    else if ( yHigh == yStart )
      xStart = xHigh;
    else {
      SFindRootIn in;
      FreqDiffParamStruc par;
      in.xmax = xHigh;
      in.xmin = xLow;
      in.xacc = ACCURACY;
      in.function = FreqDiff;
      par.c = c;
      par.b = b;
      par.y0 = yStart;
      TRY( LALSBisectionFindRoot( stat->statusPtr, &xStart, &in,
				  (void *)( &par ) ), stat );
    }
  }

  /* Compute initial dimensionless time, and record actual initial
     frequency in case it is different. */
  t = pow( xStart, -8.0 );
  FREQ( yStart, xStart );
  if ( yStart >= yMax ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EFBAD,
	   GENERATEPPNINSPIRALH_MSGEFBAD );
  }
  params->fStart = yStart*fFac;

  /*******************************************************************
   * GENERATE WAVEFORM                                               *
   *******************************************************************/

  /* Set up data pointers and storage. */
  if ( buffer ) {
    here = head =
      (PPNInspiralBuffer *)LALMalloc( sizeof(PPNInspiralBuffer) );
    if ( !here ) {
      ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
	     GENERATEPPNINSPIRALH_MSGEMEM );
    }
    here->next = NULL;
    a = here->a;
    phi = here->phi;
    nMax = (UINT4)( -1 );
    if ( params->lengthIn > 0 )
      nMax = params->lengthIn;
    nNext = BUFFSIZE;
    if ( nNext > nMax )
      nNext = nMax;
  } else {
    a = output->a->data->data;
    phi = output->phi->data->data;
    nMax = output->a->data->length;
    if ( ( params->lengthIn > 0 ) && ( nMax > params->lengthIn ) )
      nMax = params->lengthIn;
    nNext = nMax;
  }

  /* Start integrating!  Inner loop exits each time a new buffer is
     required.  Outer loop has no explicit test; when a termination
     condition is met, we jump directly from the inner loop using a
     goto statement.  All goto statements jump to the terminate: label
     at the end of the outer loop. */
  n = 0;
  dyMax = 0.0;
  y = yOld = 0.0;
  x = xStart;
  while ( 1 ) {
    while ( n < nNext ) {
      REAL4 f2a;
      REAL4 phase = 0.0;

      /* Check if we're still in a valid PN regime. */
      if ( x > xMax ) {
	params->termCode = GENERATEPPNINSPIRALH_EPNFAIL;
	params->termDescription = GENERATEPPNINSPIRALH_MSGEPNFAIL;
	goto terminate;
      }

      /* Compute the frequency.  Note that this also computes the
         global variables x2 and x3, which may be used later. */
      FREQ( y, x );
      if ( y > yMax ) {
	params->termCode = GENERATEPPNINSPIRALH_EFSTOP;
	params->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
	goto terminate;
      }
      if ( y < yOld ) {
	params->termCode = GENERATEPPNINSPIRALH_EFNOTMON;
	params->termDescription = GENERATEPPNINSPIRALH_MSGEFNOTMON;
	goto terminate;
      }
      if ( y - yOld > dyMax )
	dyMax = y - yOld;

      /* Compute the amplitude from the frequency. */
      f2a = pow( f2aFac*y, TWOTHIRDS );
      *(a++) = apFac*f2a;
      *(a++) = acFac*f2a;

      /* Compute the phase. */
      if ( b0 )
	phase += d0;
      if ( b1 )
	phase += d1*x;
      if ( b2 )
	phase += d2*x2;
      if ( b3 )
	phase += d3*x3;
      if ( b4 )
	phase += d4*x3*x;
      if ( b5 )
	phase += d5*log(t)*x3*x2;
      phase *= t*x3*etaInv;
      *(phi++) = phiC - phase;

      /* Increment the timestep. */
      n++;
      t += dt;
      yOld = y;
      if ( t < 0.0625 ) {
	params->termCode = GENERATEPPNINSPIRALH_ERTOOSMALL;
	params->termDescription = GENERATEPPNINSPIRALH_MSGERTOOSMALL;
	goto terminate;
      }
      x = pow( t, -0.125 );
    }

    /* We've either filled the buffer (if we were buffering) or we've
       exceeded the maximum length.  If the latter, we're done! */
    if ( n >= nMax ) {
      params->termCode = GENERATEPPNINSPIRALH_ELENGTH;
      params->termDescription = GENERATEPPNINSPIRALH_MSGELENGTH;
      goto terminate;
    }

    /* Otherwise, allocate the next buffer. */
    here->next =
      (PPNInspiralBuffer *)LALMalloc( sizeof(PPNInspiralBuffer) );
    here = here->next;
    if ( !here ) {
      FREELIST( head );
      ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
	     GENERATEPPNINSPIRALH_MSGEMEM );
    }
    here->next = NULL;
    a = here->a;
    phi = here->phi;
    nNext += BUFFSIZE;
    if ( nNext > nMax )
      nNext = nMax;
  }

  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  /* The above loop only exits by triggering one of the termination
     conditions, which jumps to the following point for cleanup and
     return. */
 terminate:

  /* First, set remaining output parameter fields. */
  params->dfdt = dyMax*fFac*output->a->deltaT;
  params->fStop = yOld*fFac;
  params->length = n;

  /* If data was being buffered, we need to allocate the output
     REAL4VectorSequence structures and pack the data into them. */
  if ( buffer ) {
    CreateVectorSequenceIn in;
    in.length = n;

    /* First, though, we should destroy anything already there. */
    if ( output->a->data ) {
      TRY( LALSDestroyVectorSequence( stat->statusPtr,
				      &( output->a->data ) ), stat );
    }
    if ( output->phi->data ) {
      TRY( LALSDestroyVectorSequence( stat->statusPtr,
				      &( output->phi->data ) ), stat );
    }

    /* Now try creating them. */
    in.vectorLength = 2;
    LALSCreateVectorSequence( stat->statusPtr, &( output->a->data ),
			      &in );
    BEGINFAIL( stat )
      FREELIST( head );
    ENDFAIL( stat );
    in.vectorLength = 1;
    LALSCreateVectorSequence( stat->statusPtr, &( output->phi->data ),
			      &in );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVectorSequence( stat->statusPtr,
				      &( output->a->data ) ), stat );
      FREELIST( head );
    } ENDFAIL( stat );

    /* Structures have been successfully allocated; now fill them.  We
       deallocate the list as we go along. */
    a = output->a->data->data;
    phi = output->phi->data->data;
    here = head;
    while ( here && ( n > 0 ) ) {
      PPNInspiralBuffer *last = here;
      UINT4 nCopy = BUFFSIZE;
      if ( nCopy > n )
	nCopy = n;
      memcpy( a, here->a, 2*nCopy*sizeof(REAL4) );
      memcpy( phi, here->phi, nCopy*sizeof(REAL4) );
      a += 2*nCopy;
      phi += nCopy;
      n -= nCopy;
      here = here->next;
      LALFree( last );
    }

    /* This shouldn't happen, but free any extra buffers in the
       list. */
    FREELIST( here );
  }

  /* So that's what to do if the data was buffered.  Otherwise, we've
     been filling the output arrays all along, so we just need to pad
     them out with zeros. */
  else if ( n < nMax ) {
    memset( a, 0, 2*( nMax - n )*sizeof(REAL4) );
    memset( phi, 0, ( nMax - n )*sizeof(REAL4) );
  }

  /* Everything's been stored and cleaned up, so there's nothing left
     to do but quit! */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
