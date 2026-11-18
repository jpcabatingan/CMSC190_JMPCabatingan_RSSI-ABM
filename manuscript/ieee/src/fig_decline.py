"""Fig: chemical load vs efficacy decline, window A. Black and white. Data: raw batch CSVs."""
import sys, csv, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R, OUT = sys.argv[1], sys.argv[2]
plt.rcParams.update({"font.family": "Liberation Serif", "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.6})
rows = lambda p: [r for r in csv.reader(open(f"{R}/{p}")) if r]
A, D, S, O = rows("A/main/batch_main.csv"), rows("A/main/batch_decline.csv"), rows("A/main/batch_saturating.csv"), rows("A/oat/batch_oat.csv")
load = lambda rs: [float(r[7]) for r in rs]
order = ["0", "160", "80", "40", "20", "10"]
pts = {n: load([r for r in A if r[2] == "chemical" and r[3] == n] + [r for r in D if r[3] == n]) for n in order}
bh = load([r for r in A if r[2] in ("biological", "hybrid")] + [r for r in S if r[2] in ("biological", "hybrid")]
          + [r for r in O if r[1] == "0" and r[2] in ("biological", "hybrid")])
none = st.mean(load([r for r in A if r[2] == "none"]))
fig, ax = plt.subplots(figsize=(3.4, 2.25))
ax.axhspan(min(bh), max(bh), color="0.82", lw=0)
ax.text(5.1, min(bh) / 1.12, "biological and hybrid, range of 14 runs", ha="right", va="top", fontsize=7)
ax.axhline(none, color="black", ls="--", lw=0.7)
ax.text(5.1, none * 1.15, "no control (mean of 2 runs)", ha="right", va="bottom", fontsize=7)
xs = range(len(order))
ax.plot(xs, [st.mean(pts[n]) for n in order], color="black", lw=0.8)
for i, n in enumerate(order):
    ax.scatter([i] * len(pts[n]), pts[n], s=10, facecolor="white", edgecolor="black", lw=0.7, zorder=3)
ax.set_yscale("log"); ax.set_ylim(1.5e4, 4e7)
ax.set_xticks(list(xs)); ax.set_xticklabels(["none", "160", "80", "40", "20", "10"])
ax.set_xlabel("Chemical efficacy halves every ... applications per farm")
ax.set_ylabel("Cumulative pest load (K-days)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(pad=0.2)
fig.savefig(f"{OUT}/fig-decline.pdf"); fig.savefig(f"{OUT}/fig-decline.png", dpi=300)
print({n: round(st.mean(pts[n])) for n in order}, round(min(bh)), round(max(bh)), round(none))
