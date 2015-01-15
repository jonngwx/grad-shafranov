.PHONY: all

all: Final_Report

Final_Report: Final_Report.tex
	pdflatex $^
	rm $@.log $@.aux

.PHONY: clean

clean: 
	$(RM) Final_Report.pdf
