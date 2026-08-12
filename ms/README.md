# Main manuscript

This directory contains the generic arXiv-style LaTeX version of
“Dynamic consequences of declining vaccination coverage in a heterogeneous
world.” The initial source was converted from the project Google document:

<https://docs.google.com/document/d/1O6cOvPdR1GhHBe-c8_iUdeBdXx0IcBT6NO5Gun0MfDU/edit?tab=t.0>

The conversion preserves the manuscript's wording and order. Google Docs
comments are retained as attributed author-comment macros (`\david`,
`\emily`, `\jess`, and `\bryan`). David's concrete textual suggestions are
shown as blue strikeouts (`\stkout`) and red replacement text (`\dsugg`) so
that proposed edits remain distinct from discussion comments. The default
build is deliberately journal-neutral and includes light-grey line numbers
for review.

The references are maintained in `measles_honeymoon.bib` and processed with
BibTeX. The local `measles_honeymoon.bst` is the Vancouver bibliography style,
included locally to make the arXiv build self-contained. Complete author
metadata remains in the BibTeX database; the style abbreviates long displayed
author lists according to Vancouver conventions. DOIs are printed for
publications that have them, and URLs are printed for
the remaining online sources; every bibliography entry contains at least one
of these links. From this directory:

- `make` builds `measles_honeymoon.pdf`.
- `make measles_honeymoon.open` builds and opens the PDF on macOS.
- `make measles_honeymoon.fresh` removes this manuscript's generated files
  and PDF without touching any other PDF.
