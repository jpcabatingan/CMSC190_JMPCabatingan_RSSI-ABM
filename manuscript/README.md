# Manuscript

| Folder | Format | Contents |
|---|---|---|
| `sp1-full/` | ICS thesis class (`icsthesis.cls`) | Full manuscript: front matter, Chapters 1 to 4, literature cited |
| `ieee/` | IEEE conference (`IEEEtran`) | Journal paper |

Compiled PDFs are in each folder (`main.pdf`) and in `../presentation/`.

## Build

```
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

`sp1-full/` uses `apacite`, `lipsum` and `longtable`. `ieee/` uses `IEEEtran` and `xurl`.

## Figures and numbers

`ieee/src/` holds the scripts behind the journal paper's results. Each takes the results folder as its first argument:

```
python ieee/src/recompute.py ../code/rssi_prototype/results
python ieee/src/hindcast.py ../code/rssi_prototype/results
python ieee/src/fig_decline.py ../code/rssi_prototype/results ieee/figures/fig-decline.pdf
```
