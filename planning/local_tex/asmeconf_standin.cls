% Stand-in for asmeconf v1.47, for a local layout check on TeX Live 2018 only. Same page geometry, column
% separation, 10 pt Times and small bold sans captions. Not the real class.
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{asmeconf}[2026/11/27 local stand-in]
\DeclareOption*{}
\ProcessOptions\relax
\LoadClass[twoside,10pt,twocolumn]{article}
\RequirePackage[letterpaper,left=0.5in,right=0.5in,top=0.5in,bottom=1.in,footskip=0.3in,columnsep=0.38in]{geometry}
\RequirePackage[T1]{fontenc}
\RequirePackage{mathptmx,helvet,textcomp}
\RequirePackage{mathtools,amssymb,bm,booktabs,graphicx,array,enumitem}
\RequirePackage[font={sf,bf,small},labelsep=colon]{caption}
\RequirePackage[numbers,sort&compress]{natbib}
\RequirePackage[hidelinks]{hyperref}
\RequirePackage{lineno}\linenumbers
\setlist{nosep,leftmargin=*}
\renewcommand\figurename{FIGURE}\renewcommand\tablename{TABLE}
\newcommand\ConfName[1]{\gdef\@confname{#1}}
\newcommand\ConfAcronym[1]{}\newcommand\ConfDate[1]{}\newcommand\ConfCity[1]{}
\newcommand\PaperNo[1]{\gdef\@paperno{#1}}
\newcommand\SetAuthors[1]{\gdef\@authors{#1}}
\newcommand\affil[1]{\textsuperscript{#1}}
\newcommand\CorrespondingAuthor[1]{\textsuperscript{,*}}
\newcommand\SetAffiliation[2]{\gdef\@affil{\textsuperscript{#1}#2}}
\newcommand\versionfootnote[1]{}
\newcommand\keywords[1]{\gdef\@keywords{#1}}
\renewcommand\maketitle{\twocolumn[{\centering\sffamily
  {\small\@confname\par}\vspace{6pt}{\bfseries\large\@paperno\par}\vspace{14pt}
  {\LARGE\bfseries\MakeUppercase{\@title}\par}\vspace{14pt}
  {\bfseries\@authors\par}\vspace{4pt}{\small\@affil\par}\vspace{20pt}}]}
\renewenvironment{abstract}{\section*{ABSTRACT}}{\par\medskip\noindent\textbf{Keywords: }film cooling, laidback fan-shaped hole, 7-7-7 database, direct numerical simulation, large-eddy simulation, RANS, turbulent Prandtl number\par}
\newenvironment{nomenclature}[1][2em]{\section*{NOMENCLATURE}\begin{list}{}{\setlength\itemsep{0pt}\setlength\parsep{0pt}\setlength\topsep{0pt}\setlength\labelwidth{#1}\setlength\leftmargin{\dimexpr#1+1em}\setlength\labelsep{1em}}}{\end{list}}
\newcommand\entry[2]{\item[#1\hfill]#2}
\newcommand\EntryHeading[1]{\item[]\hspace*{-3.6em}\textit{#1}}
\renewcommand\section{\@startsection{section}{1}{\z@}{2ex plus 1ex}{1ex}{\sffamily\bfseries\large\MakeUppercase}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}{1.5ex plus 1ex}{0.5ex}{\sffamily\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}{1.2ex plus 1ex}{0.3ex}{\sffamily\bfseries\small}}
\endinput
