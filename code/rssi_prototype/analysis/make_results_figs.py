"""Results figures for the 13 Nov adviser deck, in the deck palette. Data: results/series/*.csv."""
import csv, collections, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import spearmanr

sys.path.insert(0, "/mnt/user-data/uploads/SP/sp-tracker/figures")
from palette import NAVY, BLUE_STRONG, BLUE_PERI, SUB, LINE, ORANGE

S = "/mnt/user-data/uploads/SP/sp-tracker/gama/rssi_prototype/results/series/"
for f in font_manager.findSystemFonts():
    if "Poppins" in f:
        font_manager.fontManager.addfont(f)
plt.rcParams.update({"font.family": "Poppins", "font.size": 15, "axes.edgecolor": SUB,
                     "axes.labelcolor": NAVY, "xtick.color": SUB, "ytick.color": SUB,
                     "axes.spines.top": False, "axes.spines.right": False})

# Time scale. Keep only the two runs that reached day 1,092; days 0-266 also hold a third,
# incomplete run, so the mean is taken per run index within each day.
rows = [r for r in csv.reader(open(S + "series_hybrid_baseline.csv")) if r[0] == "hybrid_baseline"]
byday = collections.defaultdict(list)
for r in rows:
    byday[int(r[5])].append(float(r[7]))
days = sorted(byday)
mean = [sum(byday[d][:2]) / 2 for d in days]
lo = [min(byday[d][:2]) for d in days]
hi = [max(byday[d][:2]) for d in days]
sra_share = 72103.42 / 239215.62
sra_on_map = 72103.42 / 356869  # validated ha over cane cells in the input grid
cross = lambda th: next(d for d, m in zip(days, mean) if m >= th)

fig, ax = plt.subplots(figsize=(12, 6.2), dpi=200)
ax.fill_between(days, lo, hi, color=BLUE_PERI, alpha=0.35, lw=0)
ax.plot(days, mean, color=BLUE_STRONG, lw=3)
for th, lab, ls, dy, tx, ha in [(sra_share, f"SRA, 2 Oct 2026: 30.1% of planted cane (day {cross(sra_share)})", "--", 0.012, 15, "left"),
                    (sra_on_map, f"Same hectares on the cane map: 20.2% (day {cross(sra_on_map)})", ":", -0.032, 15, "left")]:
    ax.axhline(th, color=ORANGE, lw=2, ls=ls)
    ax.text(tx, th + dy, lab, ha=ha, color=ORANGE, fontsize=13, fontweight=700)
ax.set_xlabel("Model day (day 0 = 1 October)")
ax.set_ylabel("Share of cane ever infested")
ax.set_xlim(0, 1100); ax.set_ylim(0, 0.6)
ax.grid(axis="y", color=LINE, lw=1)
fig.tight_layout()
fig.savefig("r-timescale.png", transparent=True)
fig.savefig("r-timescale-white.png", facecolor="white")
print("cross", cross(sra_share), cross(sra_on_map), "end", mean[-1], days[-1])

# Per-LGU comparison
lg = collections.defaultdict(lambda: {"cane": 0, "ever": [], "sra": -1, "name": ""})
for r in csv.reader(open(S + "lgu_hybrid_baseline.csv")):
    g = lg[int(r[2])]
    g["name"], g["cane"], g["sra"] = r[3], int(r[4]), float(r[6])
    g["ever"].append(int(r[5]))
pts = [(g["sra"] / g["cane"], sum(g["ever"]) / len(g["ever"]) / g["cane"], g["name"])
       for g in lg.values() if g["cane"] > 0 and g["sra"] >= 0]
x, y, names = zip(*pts)
rho = spearmanr(x, y).correlation
print("n", len(pts), "rho", rho)

fig, ax = plt.subplots(figsize=(9, 6.75), dpi=200)
ax.scatter(x, y, s=60, color=BLUE_STRONG, alpha=0.85, lw=0)
ax.set_xlabel("SRA validated area / mapped cane area, per municipality")
ax.set_ylabel("Model: share of cane ever infested")
ax.set_title(f"{len(pts)} municipalities, Spearman {rho:.2f}", color=NAVY, fontsize=17, fontweight=700, loc="left")
ax.grid(color=LINE, lw=1)
for xi, yi, n in pts:
    if xi > 1:
        ax.annotate(n.split(" (")[0] + "\n(SRA area > mapped cane)", (xi, yi), xytext=(-10, 18),
                    textcoords="offset points", ha="right", fontsize=12, color=SUB)
fig.tight_layout()
fig.savefig("r-lgu.png", transparent=True)
fig.savefig("r-lgu-white.png", facecolor="white")
