
all:: elegant.ps elegant.pdf overview.ps overview.pdf Pelegant.ps Pelegant.pdf html

elegant.ps:
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc dictionary.tex
	elegant print_dict.ele
	latex elegant
	latex elegant
	dvips elegant -o elegant.ps

elegant.pdf: elegant.ps
	ps2pdf13 elegant.ps elegant.pdf

overview.ps: overview.tex
	rm -f overview.aux overview.dvi overview.log overview.toc 
	latex overview
	latex overview
	dvips overview -o overview.ps

overview.pdf: overview.ps
	ps2pdf13 overview.ps overview.pdf

Pelegant.ps: Pelegant.tex
	rm -f Pelegant.aux Pelegant.dvi Pelegant.log Pelegant.toc 
	latex Pelegant
	latex Pelegant
	dvips Pelegant -o Pelegant.ps

Pelegant.pdf: Pelegant.ps
	ps2pdf13 Pelegant.ps Pelegant.pdf

html:
	latex2html -local_icons elegant
	latex2html -local_icons Pelegant

clean::
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc elegant.ps elegant.pdf dictionary.tex
	rm -rf elegant
	rm -f overview.aux overview.dvi overview.log overview.toc overview.ps overview.pdf
	rm -rf overview
	rm -f Pelegant.aux Pelegant.dvi Pelegant.log Pelegant.toc Pelegant.ps Pelegant.pdf
	rm -rf Pelegant
