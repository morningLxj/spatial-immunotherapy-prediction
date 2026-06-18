.PHONY: check manuscript supplement all

check:
	Rscript code/check_source_data.R

manuscript:
	cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main_manuscript.tex

supplement:
	cd supplementary && latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error supplementary_material.tex

all: check manuscript supplement

