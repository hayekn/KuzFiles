\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{hayektex}[2024/04/25]

%OPTIONS{
\DeclareOption{center}{\newcommand{\@wantscenter}{}}
\DeclareOption{copyright}{\newcommand{\@copyright}{}}
\DeclareOption{authhead}{\newcommand{\@authhead}{}}
\DeclareOption{toc}{\newcommand{\@toc}{}}
\DeclareOption{title}{\newcommand{\@titpage}{}\newcommand{\@stitle}{0}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
%}

\DeclareRobustCommand{\smallcaps}[1]{\scshape{#1}}
\LoadClass[letterpaper,11pt]{article}

%PACKAGES{
\RequirePackage{geometry}
\RequirePackage{ifthen}
\ifthenelse{\isundefined{\@wantscenter}}{\geometry{verbose,letterpaper,top=1.45in,bottom=1.35in,left=1.15in,right=2in,asymmetric,twoside,marginpar=1.6in}}{\geometry{verbose,letterpaper,top=1.45in,bottom=1.35in,left=1.25in,right=1.25in,twoside}}
\RequirePackage[widermath,narrowiints,fullsumlimits,fullintlimits]{kpfonts}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage{tikz}
\RequirePackage{pgfplots}
\RequirePackage{graphicx}
\RequirePackage{empheq}
\RequirePackage[framemethod=TikZ, xcolor]{mdframed}
\RequirePackage{fancyhdr}
\RequirePackage{mathdots}
\RequirePackage{microtype}
\RequirePackage{stmaryrd}
\RequirePackage{titlesec}
\RequirePackage[titles]{tocloft}
\RequirePackage{textcase}
\RequirePackage{setspace}
\RequirePackage{xfrac}
\RequirePackage{mathtools}
\RequirePackage{faktor}
\RequirePackage{cancel}
\RequirePackage{mparhack}
\RequirePackage{booktabs}
\RequirePackage{multirow}
\RequirePackage{csquotes}
\RequirePackage{epigraph}
\RequirePackage{fontawesome}
\DeclareMathAlphabet{\bignum}{U}{BOONDOX-ds}{m}{n}
\newcommand{\1}{\bignum{1}}
\newcommand{\0}{\bignum{0}}
\RequirePackage{amsfonts}
\RequirePackage{amsthm}
\RequirePackage{amsxtra}
\RequirePackage{amssymb}
\RequirePackage{amsmath}
\RequirePackage{mathrsfs}
\RequirePackage[colorlinks,linkcolor=black,urlcolor=black,citecolor=black]{hyperref}
%}

%SECTION FORMATTING{
\renewcommand{\@author}{}
\newcommand{\@email}{}
\newcommand{\email}[1]{\renewcommand{\@email}{#1}}
\newcommand{\subtitle}[1]{\newcommand{\@subtitle}{#1}}
\renewcommand{\author}[1]{\renewcommand{\@author}{#1}}
\newcommand{\stitle}[1]{\renewcommand{\@stitle}{#1}}
\newcommand{\note}[1]{\newcommand{\@note}{#1}}
\titleformat{\section}{\relax\centering\large\lsstyle}{}{0.5em}{\otherscshape\MakeLowercase}
\titleformat{\subsection}{\relax\centering}{}{1em}{\large\itshape}
\titleformat{\subsubsection}{\relax}{}{0em}{\normalsize\textbf}
\titleformat{\paragraph}[runin]{}{}{}{\normalsize\smallcaps}
\newcommand{\@curpart}{}
\let\oldpart\part
\renewcommand{\part}[1]{\oldpart{#1}\renewcommand{\@curpart}{#1}}

\renewcommand{\cftsecfont}{}
\setcounter{secnumdepth}{-1}
\setcounter{tocdepth}{3}
%}

%TEXT MACROS{
\renewcommand{\\}{\par\vspace{1.2ex}\noindent}
\newcommand{\runin}[1]{\textls[50]{\otherscshape #1}}
\newcommand{\say}[1]{\margintext{\runin{#1}}}
\newcommand{\margintext}[1]{\marginpar{\setstretch{.8}\upshape{\footnotesize{#1}}}}
%}


%TITLE PAGE AND TOC{
\renewcommand{\maketitle}{
%\newgeometry{left=2in,right=2in}
\begin{center}
	\rule{10cm}{0.4pt}\vspace{1cm}\\
	\large{\smallcaps{\textls[150]{\MakeUppercase{\@title}}}}
	
	\ifthenelse{\isundefined{\@subtitle}}{}{{\vskip 1.5ex\large\smallcaps{\textls[100]{\protect\MakeLowercase{\protect\otherscshape{\@subtitle}}}}}}
	
	\ifthenelse{\isundefined{\@author}}{}{{\vskip 1.5ex\large\smallcaps{\textls[100]{\protect\MakeLowercase{\protect\otherscshape{\@author}}}}}}
	\ifthenelse{\isundefined{\@note}}{}{{\vskip 2.5ex\textit{\@note}}}
	\\\vspace{2.5ex}\rule{3.5cm}{0.4pt}
	\newpage
\end{center}
	\restoregeometry
}


\AtBeginDocument{
\pagenumbering{gobble}
\ifthenelse{\isundefined{\@titpage}}{}{\maketitle}
\ifthenelse{\isundefined{\@toc}}{}{\newgeometry{left=2in,right=2in}\tableofcontents
	\ifthenelse{\isundefined{\@copyright}}{}{\vspace{\fill}\\\textcopyright \; \@author\; \the\year{} \;\; \href{mailto:\@email}{\faEnvelopeO}}\restoregeometry\newpage
}
\pagenumbering{arabic}
\pagestyle{fancy}
\\
}

%}


%HEADERS{
\cfoot{}
\fancyhead[LE]{\smallcaps\@title}
\fancyhead[RE]{\smallcaps\thepage}
\ifthenelse{\isundefined{\@authhead}}{
	\fancyhead[RO]{\smallcaps\protect\MakeLowercase{\protect\textls[70]{\@curpart}}}
}{
	\fancyhead[RO]{\smallcaps\protect\MakeLowercase{\protect\textls[70]{\@author}}}
}
\fancyhead[LO]{\smallcaps\thepage}
%}

%ENVIORNMENTS{
\newcommand{\@theoremcolor}{Maroon}
\newcommand{\theoremcolor}[1]{\renewcommand{\@theoremcolor}{#1}}
\mdfdefinestyle{theorem}{
	roundcorner=10pt,
	linewidth=2pt,
	backgroundcolor=\@theoremcolor!10,
	linecolor=\@theoremcolor}
\mdfdefinestyle{proof}{
	roundcorner=5pt,
	linewidth=0pt,
	backgroundcolor=gray!10}

\newenvironment{theorem}{\begin{mdframed}[style=theorem]}{\end{mdframed}}
\renewenvironment{proof}[1][]{\say{\vspace{.4cm}\\Proof#1.}\begin{mdframed}[style=proof]}{\qed\end{mdframed}}
\newenvironment{proofn}[1][]{\say{\vspace{.4cm}\\Proof#1.}\begin{mdframed}[style=proof]}{\end{mdframed}}

%}



