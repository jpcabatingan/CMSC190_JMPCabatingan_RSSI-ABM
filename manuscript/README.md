# Manuscript

Two versions of the manuscript. The PDFs are also in `../presentation/`.

| Folder | Format | Contents | PDF |
|---|---|---|---|
| `sp1-full/` | ICS thesis class (`icsthesis.cls`, apacite) | Full draft: front matter, Chapters 1 to 3, draft Chapter 4 (results), literature cited | `CMSC190_JMPCabatingan_manuscript.pdf` |
| `ieee/` | IEEE conference (`IEEEtran`) | Five-page paper with results | `CMSC190_JMPCabatingan_journal.pdf` |

## Build

```
cd sp1-full   # or ieee
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

`sp1-full/` needs `apacite` and `lipsum`. `ieee/` needs `IEEEtran` and `xurl`.

## Recomputing the IEEE numbers

The scripts in `ieee/src/` read the raw batch CSVs. Pass the results folder as the first argument:

```
python ieee/src/recompute.py ../code/rssi_prototype/results
python ieee/src/hindcast.py ../code/rssi_prototype/results
python ieee/src/fig_decline.py ../code/rssi_prototype/results ieee/figures/fig-decline.pdf
```

## Status

The acknowledgement in `sp1-full/` is still a placeholder. Chapter 4 (results and discussion) is a draft for review.
