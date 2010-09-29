\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ligodcc}
\RequirePackage{array}
\RequirePackage{calc}
\RequirePackage{times}
\RequirePackage{vmargin}
\RequirePackage{color}
\newcommand\@ligodccptsize{2}
\DeclareOption{10pt}{\renewcommand\@ligodccptsize{0}}
\DeclareOption{11pt}{\renewcommand\@ligodccptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ligodccptsize{2}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass[letterpaper,1\@ligodccptsize pt]{article}

\setpapersize{USletter}
\setmargins{1in}{0.5in}{6.5in}{9in}{10pt}{26pt}{10pt}{36pt}

\def\ps@headings{%
  \def\@evenfoot{\vbox to 0in{
     {\@ligodocdraft} 
     \hfil\mbox{\small page \thepage{}} \hfil 
     {\@ligodocdraft}
     }}
  \def\@oddfoot{\vbox to 0in{
     {\@ligodocdraft} 
     \hfil\mbox{\small page \thepage{} of \pageref{ligodcc-lastpage}}\hfil 
     {\@ligodocdraft}
     }}
  \def\@oddhead{\vbox to 0in{
     \vskip -0.2in 
     \hfil{\@ligodocdraft} \hfil
     \vskip  0.1in
     \hfil\mbox{\small\@ligodocno}\hfil
     }}
  \def\@evenhead{\vbox to 0in{
     \vskip -0.2in 
     {\small \@ligodocdraft} \hfill {\small \@ligodocdraft}
     \vskip  0.1in
     \hfil\mbox{\small\@ligodocno}\hfil
     }}
}
\def\ps@titleheadings{%
  \let\@evenhead\@empty\let\@oddhead\@empty
  \def\@oddfoot{\hfil
    \mbox{\footnotesize\@rcsid Processed with \LaTeX\ on \today}\hfil}%
  \def\@evenfoot{\hfil
    \mbox{\footnotesize\@rcsid Processed with \LaTeX\ on \today}\hfil}%
}


\def\@ligodocno{\relax}
\def\@ligodoctype{\relax}
\def\@ligodocdist{\relax}
\def\@ligodocdraft{\relax}
\newcommand{\ligodoc}[2]{%
  \def\@ligodocno{LIGO-#1}%
  \def\@ligodoctype{-#2}%
}
\newcommand{\ligodocdist}[1]{\def\@ligodocdist{#1}}
\newcommand{\ligodocdraft}[1]{\def\@ligodocdraft{#1}}
\newcommand{\PreserveBackslash}[1]{\let\temp=\\#1\let\\=\temp}

\def\@rcsid{\relax}
\def\rcsid#1{\def\next##1#1{\def\@rcsid{\mbox{RCS ##1 --- }}}\next}

\renewcommand{\today}%
  {\number\year/\ifnum\month<10 0\fi\number\month/\ifnum\day<10 0\fi\number\day}

\renewcommand\maketitle{%
  \thispagestyle{titleheadings}%
  \begin{center}
    \large
    \mbox{LASER INTERFEROMETER GRAVITATIONAL WAVE OBSERVATORY}\\
    - LIGO -\\
    \normalsize
    CALIFORNIA INSTITUTE OF TECHNOLOGY\\
    MASSACHUSETTS INSTITUTE OF TECHNOLOGY\\[0.75in]
    \begin{tabular*}{5in}%
      {|*{3}{>{\PreserveBackslash\centering\hspace{0pt}}%
             p{5in/3-\tabcolsep*2-\arrayrulewidth*2/3}}|}
      \hline
      \textbf{Document Type} &
      \textbf{\@ligodocno\@ligodoctype} &
      \@date\rule[-2ex]{0pt}{5ex} \\
      \hline
      \multicolumn{3}{|p{5in-\tabcolsep*2-\arrayrulewidth*2}|}%
        {\rule{0pt}{5ex}\centering\parbox[t]{4in}%
          {\centering\large\bfseries\@title\rule[-2ex]{0pt}{0pt}}} \\
      \hline
      \multicolumn{3}{|p{5in-\tabcolsep*2-\arrayrulewidth*2}|}%
        {\rule{0pt}{3ex}\centering\@author\rule[-2ex]{0pt}{0ex}} \\
      \hline
    \end{tabular*} \\[0.5in]
    \textit{\small Distribution of this draft:}\vspace{1ex}\\
    {\color{red} \Large \@ligodocdraft} \\
    \@ligodocdist
    \vfill
    \begin{tabular}{c@{\qquad}c}
      \textbf{California Institute of Technology} &
      \textbf{Massachusetts Institute of Technology} \\
      \textbf{LIGO Project - MS 51-33} &
      \textbf{LIGO Project - MS 20B-145} \\
      \textbf{Pasadena CA 91125} &
      \textbf{Cambridge, MA 01239} \\
      Phone (626) 395-2129 & Phone (617) 253-4824 \\
      Fax (626) 304-9834 & Fax (617) 253-7014 \\
      E-mail: info@ligo.caltech.edu & E-mail: info@ligo.mit.edu
    \end{tabular}\\[1ex]
    WWW: http://www.ligo.caltech.edu/\\
  \end{center}
  \newpage
}

\AtEndDocument{\label{ligodcc-lastpage}}

\pagestyle{headings}
