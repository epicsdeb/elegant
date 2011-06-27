
all:: elegant.ps elegant.pdf overview.ps overview.pdf

elegant.ps:
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc dictionary.tex
	../O.linux-x86_64/elegant print_dict.ele
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

clean::
	rm -f elegant.aux elegant.dvi elegant.log elegant.toc elegant.ps elegant.pdf dictionary.tex
	rm -rf elegant
	rm -f overview.aux overview.dvi overview.log overview.toc overview.ps overview.pdf
	rm -rf overview
