#====================================================================#
#                                                                    #
# Copyright 2009 Mikael Granvik, Jenni Virtanen, Karri Muinonen,     #
#                Teemu Laakso, Dagmara Oszkiewicz                    #
#                                                                    #
# This file is part of OpenOrb.                                      #
#                                                                    #
# OpenOrb is free software: you can redistribute it and/or modify it #
# under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or  #
# (at your option) any later version.                                #
#                                                                    #
# OpenOrb is distributed in the hope that it will be useful, but     #
# WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                           #
#                                                                    #
# You should have received a copy of the GNU General Public License  #
# along with OpenOrb. If not, see <http://www.gnu.org/licenses/>.    #
#                                                                    #
#====================================================================#
#
# Produces documentation for the OpenOrb project.
#
# Author:  MG
# Version: 2009-07-25

include ../make.config

FILE = OpenOrb_Tutorial
PSFILE = $(FILE).ps
PDFFILE = $(FILE).pdf
DVIFILE = $(FILE).dvi
BIBTEX = bibtex
LATEX = latex
PS2PDF= ps2pdf
DVIPS = dvips
GV = ggv -scale 0
XDVI = xdvi

#----------------------

.SUFFIXES: $(SUFFIXES) .tex .dvi .ps

.tex.dvi:
	$(LATEX) $<
	$(LATEX) $<
	$(BIBTEX) $(FILE)
	$(BIBTEX) $(FILE)
	$(LATEX) $<
	$(LATEX) $<

.dvi.ps:
	$(DVIPS) -o $(PSFILE) $<

all: $(PSFILE)

dvi: $(DVIFILE)
	$(XDVI) $(DVIFILE)

pdf: $(DVIFILE)
	$(DVIPS) -Ppdf -j0 -G0 -o prepdf_$(PSFILE) $<
	$(PS2PDF) -sPAPERSIZE=a4 -dMaxSubsetPct=100 \
	-dCompatibilityLevel=1.2 -dSubsetFonts=true \
	-dEmbedAllFonts=true -sOutputFile=$(PDFFILE) \
	prepdf_$(PSFILE)
	rm -f prepdf_$(PSFILE)

clean:
	rm -f $(PSFILE) prepdf_$(PSFILE) \
	*.dvi *.log *.aux *~ core *.bbl *.blg *.out \
	*.toc

$(PSFILE): $(DVIFILE)
$(PDFFILE): $(DVIFILE)
$(DVIFILE): $(TEXFILES)
