% !TEX root = diss.tex
%% LaTeX2e class for diploma theses
%% Based on wissdoc.cls by Roland Bless, 1996-2001
%%
%% Heavily modified by Robert Schittny, 2015


\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{thesisclass}
\LoadClass[BCOR=8mm,        
           DIV=calc, 
           12pt,
%           draft=true,
           titlepage,
           numbers=noenddot,
           toc=bibliography]
           {scrbook}
%% used class options:
% BCOR: binding correction to add to the inner page margin
% DIV=calc: calculate the optimal page setup
% 12pt: font size to use
% draft: shows black boxes at the end of overfull boxes => remove for final
%   version!
% titlepage: use separate pages for the title, i.e. the following content
%   appears on a new page
% numbers=noenddot: don't put dots at the end of chapter numbers, i.e. '1.3.4'
%   instead of '1.3.4.'
% toc=bibliography: create TOC entries for bibliography

%% Set up fonts and page spread
% We want the Palatino font and enable old style figures ([osf]) to get numbers that fit
% better into the text
\usepackage[osf]{mathpazo}	

% We want to consider the header as part of the page spread, but not the footer
\KOMAoptions{headinclude=true, footinclude=false}

% Recalculate the page layout, taking into account the chosen font and the headinclude
% option
\recalctypearea

%% ----------------------------------------------------------------------------
%%                                  General settings
%% ----------------------------------------------------------------------------


% For inverse search from within pdf documents
\synctex=1
\usepackage{pdfsync}

% Delete package option 'noframe' to show frames around boxes, graphics, etc.
\usepackage[noframe]{showframe}

% To create lorem ipsum placeholder text
\usepackage{lipsum}

% For absolute positioning of text anywhere on the page (used only for the title pages)
\usepackage[absolute]{textpos}

% Suppresses warning due to non-implemented font sizes
\usepackage{anyfontsize} 

% General enhancements of typesetting, results in better gray level of the pages
%\usepackage[tracking=smallcaps]{microtype}
\usepackage[letterspace=40, tracking=smallcaps]{microtype}

% we want to compile to pdf
\usepackage[pdftex]{graphicx} 
\usepackage{tikz}

% --- Colors --- 
\usepackage{color}  % colored text, also define some colors
\usepackage{xcolor}
\definecolor{olivegreen}{rgb}{0,0.7,0}
\definecolor{gray70}{rgb}{0.3,0.3,0.3}
\xdefinecolor{darkblue}{cmyk}{1, 1, 0, 0.1}
\definecolor{RoyalBlue}{cmyk}{1, 0.50, 0, 0}
\definecolor{webgreen}{rgb}{0,0.5,0}

% --- Todo notes ---
\usepackage[textsize=small,colorinlistoftodos]{todonotes}
\newcommand{\citmissing}[1]{\todo[color=green!40]{#1}}
\newcommand{\todolow}[1]{\todo[color=yellow!40]{#1}}
\newcommand{\inltodo}[1]{\textcolor{red}{#1}}

% --- Flexible quotes ---
\usepackage{csquotes}

% --- hyperlinks ---
\usepackage{url} % URL typesetting, e.g. in references
% Create hyperlinks for the pdf version, remove colors for printing
\usepackage[raiselinks=true,
            unicode=true,
			bookmarks=true,
			bookmarksopenlevel=1,
			bookmarksopen=true,
			bookmarksnumbered=true,
			hyperindex=true,
			plainpages=false,
			pdfpagelabels=true,
			colorlinks=true,						
            linkcolor=RoyalBlue,
            urlcolor=RoyalBlue,
            citecolor=webgreen,
            linktocpage=true,
		]{hyperref}

%% ----------------------------------------------------------------------------
%%                                  Page layout
%% ----------------------------------------------------------------------------

\usepackage{textcase} % provides \MakeTextLowercase{}

% Creates lowercase small caps with some extra letter spacing
\newcommand{\spacedlowsmallcaps}[1]{\textls[80]{\textsc{\MakeTextLowercase{#1}}}}

% Creates small caps with extra letter spacing
\newcommand{\spacedsmallcaps}[1]{\textls[40]{\scshape #1}}

% --- Header and footer ---
\usepackage[headsepline]{scrlayer-scrpage}
\pagestyle{scrheadings}
\clearscrheadfoot
\ohead{\spacedlowsmallcaps{\headmark}}  % header has a separating line
\ofoot[\pagemark]{\pagemark}    % footer is only the page number, no line

% --- Section numbering ---
\setcounter{secnumdepth}{3} % Numbering also for subsubsections
\setcounter{tocdepth}{2}    % Register only subsections in table of contents

% --- Define the look of sections, subsections etc. ---
\usepackage{titlesec}
\titleformat{\part}[display]
    {\large\filcenter}
    {Part \thepart\thispagestyle{empty}}
    {1pc}
    {\LARGE\spacedlowsmallcaps}

\titleformat{\section}
    {\large\scshape}        % format for the whole output
    {\thesection}           % which label to print
    {0.75em}                % separator between label and title
    {\spacedsmallcaps}      % code before the title

\titleformat{\subsection}
    {\normalsize\itshape}
    {\thesubsection}
    {0.75em}
    {}

\titleformat{\subsubsection}
    {\normalsize\itshape}
    {\thesubsubsection}
    {0.75em}
    {}

\titleformat{\paragraph}[runin]
    {\normalsize\scshape}
    {\theparagraph}
    {1em}
    {\spacedsmallcaps}
    [\!\!\!\!\!]

% --- define look of the fancy chapter headings
\definecolor{dgray}{gray}{0.3}
\setkomafont{chapter}{\normalfont\huge\rmfamily\scshape\lsstyle}
\addtokomafont{chapter}{\color{dgray}}

\makeatletter% siehe De-TeX-FAQ
\renewcommand*{\chapterformat}{%
\begingroup% damit \unitlength-Änderung lokal bleibt
\setlength{\unitlength}{1mm}%
\begin{picture}(20,30)(0,5)%
\setlength{\fboxsep}{0pt}%
\put(20,15){\line(1,0){\dimexpr
\textwidth-20\unitlength\relax\@gobble}}%
\put(6.5,0){\makebox(20,20)[r]{%
    \fontsize{28\unitlength}{1}\color{dgray}\selectfont\textls[-40]{\thechapter}
}}%
\put(20,15){\makebox(\dimexpr
\textwidth-20\unitlength\relax\@gobble,\ht\strutbox\@gobble)[l]{%
\ \normalsize\rmfamily\normalfont\color{black}\chapapp~\thechapter\autodot
}}%
\end{picture} % <- Leerzeichen ist hier beabsichtigt!
\endgroup
}%
\makeatother


\renewcommand{\descriptionlabel}[1]{\hspace*{\labelsep}\spacedlowsmallcaps{#1}}

\titlespacing*{\section}{0pt}{1.25\baselineskip}{1\baselineskip} 
\titlespacing*{\subsection}{0pt}{1.25\baselineskip}{1\baselineskip}
\titlespacing*{\paragraph}{0pt}{1\baselineskip}{1\baselineskip}

% --- Style of Table of Contents ---
\usepackage{titletoc}
\titlecontents{part}[1.5em]{\large\vspace{2.0em}}
    {\contentslabel{1.5em}\scshape\lsstyle}
    {\scshape\lsstyle\hspace*{-1.5em}}
    {\hspace{1em minus 1em}\normalsize\contentspage}
    [\vspace{0.25em}]

\titlecontents{chapter}[1.5em]{\vspace{0.45em}}
    {\contentslabel{1.5em}\scshape\lsstyle}
    {\scshape\lsstyle\hspace*{-1.5em}}
    {\hspace{1em minus 1em}\contentspage}
    [\vspace{0.25em}]

\titlecontents{section}[4em]{\vspace{0.1em}}
    {\contentslabel{2.5em}}
    {\hspace*{-2.5em}}
    {\hspace{1em minus 1em}\contentspage}
    [\vspace{0.1em}]

\titlecontents{subsection}[7.25em]{\vspace{0.1em}}
    {\contentslabel{3.25em}}
    {\hspace*{-3.25em}}
    {\hspace{1em minus 1em}\contentspage}
    [\vspace{0.1em}]

% --- option for chapters that don't belong to a part
\makeatletter
\newcommand{\ChapterOutsidePart}{%
   \def\toclevel@chapter{-1}\def\toclevel@section{0}\def\toclevel@subsection{1}}
\newcommand{\ChapterInsidePart}{%
   \def\toclevel@chapter{0}\def\toclevel@section{1}\def\toclevel@subsection{2}}
\makeatother

% --- define the look of footnotes ---
% (all lines indented by the same amount)
\deffootnote{1em}{1em}{\textsuperscript{\thefootnotemark}\ }



%% ----------------------------------------------------------------------------
%%                              Language and encoding
%% ----------------------------------------------------------------------------

\usepackage[ngerman]{babel}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\DeclareUnicodeCharacter{00A0}{~}

%% ----------------------------------------------------------------------------
%%                                  Math typesetting
%% ----------------------------------------------------------------------------

\usepackage[intlimits]{amsmath}
\usepackage{amssymb}
\usepackage{mathtools}
\usepackage{siunitx}
% Setup siunitx
\sisetup{%detect-mode,
    %mode = text,
    group-digits = false,
    input-signs = ,
    input-symbols = ()[]-+,
    input-open-uncertainty = ,
    input-close-uncertainty = ,
    table-align-text-post = false,
    range-phrase = --,
    range-units = single}
\usepackage{xfrac}

% typesetting of chemical formulas
\usepackage[version=4]{mhchem}

%% ----------------------------------------------------------------------------
%%                               Bibliography
%% ----------------------------------------------------------------------------

\usepackage[backend=biber, % use biber backend instead of bibtex
            style=phys, % OSA-compatible citation style
            biblabel=brackets, % put numbers in the bibliography in brackets
            backref=true % create back-references to pages in the bibliography
          ]{biblatex}

\DefineBibliographyStrings{english}{%
  backrefpage = {cited on page},
  backrefpages = {cited on pages},
}

% citation command for just the citation number without any brackets
\DeclareCiteCommand{\barecite}%
  {\usebibmacro{cite:init}%
   \usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
   \usebibmacro{cite:comp}}
  {}
  {\usebibmacro{cite:dump}%
   \usebibmacro{postnote}}

%% ----------------------------------------------------------------------------
%%                               Graphics and tables
%% ----------------------------------------------------------------------------

\usepackage{rotating}

%% Define standard dimensions for figures
\newlength{\pwidthOne} 
\setlength{\pwidthOne}{0.975\textwidth} 
\newlength{\pwidthTwo} 
\setlength{\pwidthTwo}{0.480\textwidth} 
\newlength{\pwidthThree} 
\setlength{\pwidthThree}{0.315\textwidth} 
\newlength{\pwidthThreeTwo} 
\setlength{\pwidthThreeTwo}{0.643\textwidth} 
\newlength{\pwidthFour} 
\setlength{\pwidthFour}{0.232\textwidth} 


% --- Better treatment of figures ---
% See p.105 of "TeX Unbound" for suggested values.
% See pp. 199-200 of Lamport's "LaTeX" book for details.
%   General parameters, for ALL pages:
\renewcommand{\topfraction}{0.9}	% max fraction of floats at top
\renewcommand{\bottomfraction}{0.8}	% max fraction of floats at bottom
%   Parameters for TEXT pages (not float pages):
\setcounter{topnumber}{3}
\setcounter{bottomnumber}{3}
\setcounter{totalnumber}{4}     % 2 may work better
\setcounter{dbltopnumber}{2}    % for 2-column pages
\renewcommand{\dbltopfraction}{0.9}	% fit big float above 2-col. text
\renewcommand{\textfraction}{0.07}	% allow minimal text w. figs
%   Parameters for FLOAT pages (not text pages):
\renewcommand{\floatpagefraction}{0.85}	% require fuller float pages
% N.B.: floatpagefraction MUST be less than topfraction !!
\renewcommand{\dblfloatpagefraction}{0.7}	% require fuller float pages
\setlength{\textfloatsep}{12pt plus 2pt minus 4pt}


% --- Caption settings ---
% make figure and table labels small, no caption indentation
\usepackage[normal,small,tableposition=top]{caption}
% create a small margin left and right of the caption to help further discern it from the
% main text
\setlength{\captionmargin}{1em}

% use short forms for figure and table labeling
%\addto\captionsenglish{
%\renewcommand{\figurename}{Figure}%
%\renewcommand{\tablename}{Table}%
%}

% --- Other graphics stuff
\usepackage{placeins} % force placement of figures with \FloatBarrier
\usepackage{float} % Force placement of figures with [H] (don't use if not absolutely necessary!)

\usepackage{subcaption}
\usepackage{sidecap} % offers SCfigure environment for captions beside the figure

\DeclareGraphicsExtensions{.svg,.pdf,.png,.jpg}

% --- Table stuff

% Enhanced tabular environments
\usepackage{array}
\usepackage{booktabs}
\usepackage{multirow}

% short version of the tablenum command (from siunitx)
\newcommand{\tn}[1]{\tablenum{#1}}

%% ----------------------------------------------------------------------------
%%                                  Own macros
%% ----------------------------------------------------------------------------

% Chapter abstract after chapter title
\newenvironment{chapterabstract}{\itshape}{}

\newcommand{\changefont}[3]{\fontfamily{#1} \fontseries{#2} \fontshape{#3} \selectfont}

% Creats an uneven blank page
\newcommand{\blankpage}{
 \clearpage{\pagestyle{empty}\cleardoublepage}
}

% automatically re-export svg files with LaTeX for all text upon change
\newcommand{\executeiffilenewer}[3]{\ifnum\pdfstrcmp{\pdffilemoddate{#1}}{\pdffilemoddate{#2}}>0{\immediate\write18{#3}}\fi}
\newcommand{\includesvg}[1]{\executeiffilenewer{#1.svg}{#1.pdf}{inkscape -z -C
--file=#1.svg --export-pdf=#1.pdf --export-latex}\small{\input{#1.pdf_tex}}}
  
\newcommand{\needcite}{\textcolor{blue}{[!citation!]}}
\newcommand{\needref}{\textcolor{red}{[!reference!]}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% consistent and not too big vector and tensor arrows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I want the arrows to be from the mathabx package without completely loading it (it would
% replace many other smybols too), so we define the font manually
\DeclareFontFamily{U}{matha}{\hyphenchar\font45}
\DeclareFontShape{U}{matha}{m}{n}{
      <5> <6> <7> <8> <9> <10> gen * matha
      <10.95> matha10 <12> <14.4> <17.28> <20.74> <24.88> matha12
      }{}
\DeclareSymbolFont{matha}{U}{matha}{m}{n}

% Define a subset character from that font (from mathabx.dcl)
\DeclareMathSymbol{\varleftarrow}{3}{matha}{"D0}
\DeclareMathSymbol{\varrightarrow}{3}{matha}{"D1}
\DeclareMathSymbol{\varleftrightarrow}{3}{matha}{"D8}

%% new \oset macro for overset with reduced distance
\makeatletter
\newcommand{\oset}[3][0ex]{%
  {\mathop{#3}\limits^{
    \vbox to#1{\kern-1.5\ex@
    \hbox{$\scriptstyle#2$}\vss}}}}
\makeatother

% new vector and tensor macros
\renewcommand{\vec}[1]{\oset{\varrightarrow}{#1}}
\newcommand{\tens}[1]{\oset{\varleftrightarrow}{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% protected flexible-width half-space
\newcommand{\pfhs}{\nobreak\hspace{.16667em plus .08333em}}

\newcommand{\eg}{e.g.\@}

% Upright comma and period for equations, with a small distance to the math part
\newcommand{\mcomma}{\ \text{,}}
\newcommand{\mperiod}{\ \text{.}}

% transposition in math mode
\newcommand{\transp}{\text{T}}

% Upright d for differentials
\newcommand{\md}{\,\mathrm{d}}

% Symbols I often use
\newcommand{\mus}{\mu_\text{s}}
\newcommand{\mua}{\mu_\text{a}}
\newcommand{\muam}{\mu_\text{a,m}}
\newcommand{\muap}{\mu_\text{a,p}}
\newcommand{\muex}{\mu_\text{ex}}
\newcommand{\mut}{\mu_\text{t}}
\newcommand{\mueff}{\mu_\text{eff}}
\newcommand{\ls}{l_\text{s}}
\newcommand{\lssur}{l_\text{s,\num{0}}}
\newcommand{\lsshl}{l_\text{s,\num{2}}}
\newcommand{\lt}{l_\text{t}}
\newcommand{\Ltr}{L(\!\vec{r}\!,\hat{s},t)}
\newcommand{\parspace}{\quad}
\newcommand{\Je}{\vec{J}_\text{\!\!e}}
\newcommand{\Jp}{\vec{J}}
\newcommand{\np}{n_\text{p}}
\newcommand{\nph}{n_\text{p,hom}}
\newcommand{\npcw}{n_\text{p,cw}}
\newcommand{\nphcw}{n_\text{p,hom,cw}}
\newcommand{\ze}{z_\text{e}}
\newcommand{\zs}{z_\text{s}}
\newcommand{\rp}{\vec{r}_\text{\!p}}
\newcommand{\vp}{\vec{v}_\text{\!p}}
\newcommand{\frodo}{\textsc{frodo}}
\newcommand{\tio}{\ce{TiO2}}
\newcommand{\taudiff}{\tau_\text{diff}}
\newcommand{\accuflect}{Accuflect\textregistered\,B6}
\newcommand{\Np}{N_\text{p}}
\newcommand{\tocred}[1]{\textcolor{red}{#1}}
\newcommand{\tocorange}[1]{\textcolor{orange}{#1}}
\newcommand{\tocgreen}[1]{\textcolor{green}{#1}}
\newcommand{\aDobs}{\overline{\Delta}_\text{obs}}
\newcommand{\aDclk}{\overline{\Delta}_\text{clk}}
\newcommand{\aDclks}{\overline{\Delta}_\text{clk,2}}
