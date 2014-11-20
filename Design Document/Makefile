.PHONY: all

all: Design_Document

Design_Document: Design_Document.tex
	pdflatex $^
	rm $@.log $@.aux

.PHONY: clean

clean: 
	$(RM) Design_Document.pdf
