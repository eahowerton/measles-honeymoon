MS = measles_honeymoon

export RED  := $(shell tput setaf 1)
export BLUE := $(shell tput setaf 4)
export BOLD := $(shell tput bold)
export NC   := $(shell tput sgr0)

.PRECIOUS: %.pdf %.tex

SOURCES = $(MS).tex $(MS)_preamble.tex $(MS)_main.tex $(MS).bib $(MS).bst
FIGURES = figures/figure1_global_mcv1_coverage.png \
	../figures/honeymoon_SIR_with_linear_heatmap.pdf \
	../figures/release_examples_age_AI.pdf \
	../figures/honeymoon_age_full_start95.pdf \
	../figures/empirical_vax_declines_v3_AI.pdf

LATEXMK      := latexmk
LATEXMKFLAGS := -pdf -bibtex -silent -use-make

.PHONY: all clean veryclean superclean fresh FORCE

all: $(MS).pdf

%.open: %.pdf
	open "$<"

$(MS).pdf: $(SOURCES) $(FIGURES) FORCE
	@echo $${BOLD}$${RED}tex to pdf ...$${NC}
	$(LATEXMK) $(LATEXMKFLAGS) "$(MS)"
	@echo $${BOLD}$${BLUE}`date`$${NC}

FORCE:
	@:

clean:      $(MS).clean
veryclean:  $(MS).veryclean
superclean: $(MS).superclean
fresh:      $(MS).fresh

%.clean: FORCE
	@rm -f "$*.aux" "$*.tikz" "$*.out" "$*.log" "$*.blg" "$*.toc"

%.veryclean: %.clean FORCE
	@rm -f "$*.nav" "$*.snm" "$*.fls" "$*.fdb_latexmk" "$*-tikzDictionary" "$*.synctex.gz"

%.superclean: %.veryclean FORCE
	@rm -f "$*.bbl" "$*.bcf" "$*.run.xml" "$*.ent"

%.fresh: %.superclean FORCE
	@rm -f "$*.pdf"
	@if [ -f "$*.Rnw" ]; then rm -f "$*.tex"; fi
