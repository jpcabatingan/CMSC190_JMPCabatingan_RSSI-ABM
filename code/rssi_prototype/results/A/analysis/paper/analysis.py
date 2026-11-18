"""Analysis of the window A seeded runs (main comparison, sensitivity, wind, decline and hindcast).
Run:  python analysis.py   (from this folder). Reads ../../{main,oat,wind}/*.csv, writes numbers.md and fig-*.png here.
Needs pandas and matplotlib.
"""
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, "..", "..")
COLS = ["experiment", "oat_case", "strategy", "n_half", "seed", "share_ever", "share_end", "load",
        "apps_chem", "apps_bio", "cane_total", "harvest_count", "max_D"]
N_FARMS = 10611   # farmer agents in window A (last column of a window A weekly series, test run 13 Nov); not in these CSVs
YEARS = 3.0       # 1,095 days after seeding
ORDER = ["none", "chemical", "biological", "hybrid"]
CB = {"none": "#000000", "chemical": "#D55E00", "biological": "#009E73", "hybrid": "#0072B2"}   # Okabe-Ito
OAT = {1: "settlement survival 0.01", 2: "settlement survival 0.10", 3: "neighbour transmission x0.25",
       4: "neighbour transmission x4", 5: "long jump 0", 6: "long jump 0.2", 7: "parasitism 0", 8: "parasitism 0.668",
       9: "bio:chem ratio 0.8", 10: "bio:chem ratio 1.0", 11: "level thresholds x0.5", 12: "level thresholds x2",
       13: "trigger T = 2", 14: "trigger T = 3", 15: "re-treatment 7 d", 16: "re-treatment 28 d",
       17: "wind bias 0.5", 18: "wind bias 1.0"}
PARAMS = {"settlement survival": (1, 2), "neighbour transmission": (3, 4), "long jump": (5, 6), "parasitism": (7, 8),
          "bio:chem ratio": (9, 10), "level thresholds": (11, 12), "trigger level": (13, 14),
          "re-treatment interval": (15, 16), "wind bias": (17, 18)}

def rd(p):
    return pd.read_csv(os.path.join(A, p), header=None, names=COLS)

main, oat, wind, sat = rd("main/batch_main.csv"), rd("oat/batch_oat.csv"), rd("wind/batch_wind_direction.csv"), rd("main/batch_saturating.csv")
for d in (main, oat, wind, sat):
    d["ln"] = np.log(d.load); d["apps"] = d.apps_chem + d.apps_bio; d["cane"] = d.cane_total / d.harvest_count

def lab(s, n):
    return s if s in ("none", "biological") else f"{s} (n_half {int(n) or 'none'})"

out = ["# Window A seeded runs: numbers", "", "All values are computed from the raw CSVs by analysis.py.", ""]
def P(*a): out.extend(a)
def tbl(df, fmt=None):
    df = df.copy()
    for c in df.columns:
        if fmt and c in fmt: df[c] = df[c].map(fmt[c])
    return ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] + \
           ["| " + " | ".join(map(str, r)) + " |" for r in df.itertuples(index=False)] + [""]

# ---------- 1 Reproduce section 1 of summary.md
g = main.groupby(["strategy", "n_half"])
t = g.agg(n=("load", "size"), mean=("load", "mean"), sd=("load", "std"), lo=("load", "min"), hi=("load", "max"),
          ever=("share_ever", "mean"), end=("share_end", "mean"), chem=("apps_chem", "mean"), bio=("apps_bio", "mean"),
          cane=("cane", "mean"), maxD=("max_D", "max")).reset_index()
t["config"] = [lab(s, n) for s, n in zip(t.strategy, t.n_half)]
t = t.sort_values("mean").reset_index(drop=True)
P("## 1. Main comparison (recomputed)", "")
P(*tbl(t[["config", "n", "mean", "sd", "lo", "hi", "ever", "end", "chem", "bio", "cane", "maxD"]],
       {"mean": "{:,.0f}".format, "sd": "{:,.0f}".format, "lo": "{:,.0f}".format, "hi": "{:,.0f}".format, "ever": "{:.3f}".format,
        "end": "{:.3f}".format, "chem": "{:,.0f}".format, "bio": "{:,.0f}".format, "cane": "{:.2f}".format, "maxD": "{:.3f}".format}))
sm = os.path.join(A, "analysis", "summary.md")
if os.path.exists(sm):
    txt = open(sm, encoding="utf-8").read()
    sec = txt[txt.index("## 1."):txt.index("## 2.")]
    old = {m.group(1).strip(): float(m.group(2).replace(",", "")) for m in re.finditer(r"\| \d+ \| (.+?) \| \d+ \| ([\d,]+) \(", sec)}
    diffs = {c: abs(old[c] - round(m)) for c, m in zip(t.config, t["mean"]) if c in old}
    P(f"Check against summary.md section 1: {len(diffs)} of {len(t)} configurations matched by name; "
      f"largest difference in mean load = {max(diffs.values()):.0f} K-days (rounding only if 0 or 1).", "")

# ---------- 2 Replicate variability
P("## 2. Replicate variability", "")
v = g.agg(cv_load=("load", lambda x: x.std() / x.mean()), rr_load=("load", lambda x: (x.max() - x.min()) / x.mean()),
          rr_ever=("share_ever", lambda x: (x.max() - x.min()) / x.mean()),
          rr_apps=("apps", lambda x: (x.max() - x.min()) / x.mean() if x.mean() > 0 else np.nan)).reset_index()
v["config"] = [lab(s, n) for s, n in zip(v.strategy, v.n_half)]
P("Main comparison, n = 2 per configuration. CV = SD / mean; relative range = (max - min) / mean.", "")
P(*tbl(v[["config", "cv_load", "rr_load", "rr_ever", "rr_apps"]], {c: "{:.1%}".format for c in ["cv_load", "rr_load", "rr_ever", "rr_apps"]}))

# pooled within-group SD of ln(load), by strategy, over every replicate pair in the data (main 8, wind 8, saturating 4)
groups = [(k[0], d) for k, d in main.groupby(["strategy", "n_half"])] + [(k[1], d) for k, d in wind.groupby(["oat_case", "strategy"])] + \
         [(k, d) for k, d in sat.groupby("strategy")]
ss, df_ = {}, {}
for s, d in groups:
    ss[s] = ss.get(s, 0) + ((d.ln - d.ln.mean()) ** 2).sum(); df_[s] = df_.get(s, 0) + len(d) - 1
S = {s: float(np.sqrt(ss[s] / df_[s])) for s in ss}
P("Pooled within-pair SD of ln(load) by strategy, from all 20 replicate pairs (main 8, wind 8, saturating 4). "
  "For small values this is about the CV.", "")
P(*tbl(pd.DataFrame({"strategy": ORDER, "pairs": [df_[s] for s in ORDER], "sd_ln": [S[s] for s in ORDER]}), {"sd_ln": "{:.3f}".format}))

# baseline pool: loss does not feed back on the pest (model lines 531-535), so the saturating-damage runs and OAT case 0
# are extra baseline replicates for pest load (n_half 20 for chemical and hybrid).
b_main = main[(main.n_half == 20) | main.strategy.isin(["none", "biological"])]
base = pd.concat([b_main, sat, oat[oat.oat_case == 0]])
BM = base.groupby("strategy").load.agg(["size", "mean", "std", "min", "max"]).loc[ORDER]
BM["cv"] = BM["std"] / BM["mean"]
P("Baseline pool for pest load (main at n_half 20 + saturating-damage runs + OAT case 0; damage does not affect the pest, "
  "so these are replicates of one process):", "")
P(*tbl(BM.reset_index(), {"mean": "{:,.0f}".format, "std": "{:,.0f}".format, "min": "{:,.0f}".format, "max": "{:,.0f}".format, "cv": "{:.1%}".format}))

# adjacent pairs in the main ranking: ratio of means vs noise
P("Adjacent configurations in the main ranking: ratio of mean loads and whether it exceeds the noise band "
  "exp(2 x sqrt(s_a^2/2 + s_b^2/2)) from the pooled SDs (n = 2 each).", "")
rows = []
for i in range(len(t) - 1):
    a, b = t.loc[i], t.loc[i + 1]
    band = float(np.exp(2 * np.sqrt(S[a.strategy] ** 2 / 2 + S[b.strategy] ** 2 / 2)))
    r = b["mean"] / a["mean"]
    rows.append([a.config, b.config, r, band, "distinct" if r > band else "tie"])
P(*tbl(pd.DataFrame(rows, columns=["lower load", "higher load", "ratio", "noise band", "verdict"]), {"ratio": "{:.2f}".format, "noise band": "{:.2f}".format}))

# ---------- 3 Break-even of chemical vs biological/hybrid
P("## 3. Efficacy decline: where chemical would match biological/hybrid", "")
c0, c20, c10 = [t[(t.strategy == "chemical") & (t.n_half == n)]["mean"].iloc[0] for n in (0, 20, 10)]
allrows = pd.concat([main, sat, oat[oat.oat_case == 0]])
hyb_all = allrows[(allrows.strategy == "hybrid") & (allrows.apps_chem == 0)]
bio_all = base[base.strategy == "biological"]
bh = pd.concat([hyb_all, bio_all]).load
target = bh.mean()
P(f"Chemical mean load: no decline {c0:,.0f}; n_half 20 {c20:,.0f} ({c20 / c0:.1f}x); n_half 10 {c10:,.0f} ({c10 / c0:.1f}x).",
  f"Biological and hybrid baseline runs pooled (n = {len(bh)}): mean {target:,.0f}, range {bh.min():,.0f}-{bh.max():,.0f}.",
  f"Ratios to that pool: chemical no decline {c0 / target:.2f}; n_half 20 {c20 / target:.1f}x; n_half 10 {c10 / target:.1f}x.", "")
f_log = (np.log(target) - np.log(c0)) / (np.log(c20) - np.log(c0)); f_lin = (target - c0) / (c20 - c0)
P("Break-even between 'no decline' (decline rate 1/n_half = 0) and n_half 20 (rate 0.05), two interpolations over only two points:",
  f"- log-linear in load: rate {0.05 * f_log:.4f}, n_half about {1 / (0.05 * f_log):,.0f}",
  f"- linear in load: rate {0.05 * f_lin:.5f}, n_half about {1 / (0.05 * f_lin):,.0f}",
  "The two differ by a factor of several: the data only bound the break-even at n_half > 20.", "")

# ---------- 4 Hybrid vs biological
P("## 4. Hybrid vs biological", "")
every = pd.concat([main.assign(src="main"), sat.assign(src="saturating"), wind.assign(src="wind"), oat.assign(src="oat")])
h = every[every.strategy == "hybrid"]
nz = h[h.apps_chem > 0]
P(f"Hybrid runs in all files: {len(h)}; with any chemical application: {len(nz)} "
  f"({', '.join(f'{r.src} case {int(r.oat_case)}: {int(r.apps_chem):,} chemical, {int(r.apps_bio):,} biological' for r in nz.itertuples())}).",
  f"Max cell density in the {len(hyb_all)} baseline hybrid runs: {hyb_all.max_D.min():.3f}-{hyb_all.max_D.max():.3f} K; "
  f"farms are 2x2 cells and Level 3 needs a FARM MEAN of 0.11 K, so one cell at that density gives a farm mean near {hyb_all.max_D.max() / 4:.3f} K if its neighbours are clean.", "")
o10 = oat[oat.oat_case == 10].set_index("strategy").load
P(f"Consistency check: with bio:chem ratio 1.0 (case 10) biological is the same process as chemical without decline: "
  f"biological case 10 = {o10['biological']:,.0f} vs chemical no decline = {c0:,.0f} (ratio {o10['biological'] / c0:.2f}).", "")
hm, bm = hyb_all.load, bio_all.load
se = np.sqrt(hyb_all.ln.var() / len(hm) + bio_all.ln.var() / len(bm))
P(f"Baseline, no chemicals used: hybrid n = {len(hm)} mean {hm.mean():,.0f} (CV {hm.std() / hm.mean():.1%}); "
  f"biological n = {len(bm)} mean {bm.mean():,.0f} (CV {bm.std() / bm.mean():.1%}); ratio hybrid/biological {hm.mean() / bm.mean():.3f}; "
  f"difference of mean ln(load) = {hyb_all.ln.mean() - bio_all.ln.mean():+.3f}, SE {se:.3f}.", "")
P("OAT cases where hybrid could differ (single runs):", "")
rows = []
for c in (11, 12, 13, 14, 15, 16):
    o = oat[oat.oat_case == c].set_index("strategy")
    rows.append([c, OAT[c], int(o.apps_chem["hybrid"]), int(o.apps_bio["hybrid"]), o.load["hybrid"] / o.load["biological"],
                 o.load["hybrid"] / o.load["chemical"], o.max_D["hybrid"]])
P(*tbl(pd.DataFrame(rows, columns=["case", "assumption", "hybrid chem apps", "hybrid bio apps", "hybrid/biological", "hybrid/chemical", "hybrid max_D"]),
       {"hybrid chem apps": "{:,}".format, "hybrid bio apps": "{:,}".format, "hybrid/biological": "{:.2f}".format, "hybrid/chemical": "{:.2f}".format, "hybrid max_D": "{:.2f}".format}))

# ---------- 5 OAT screen
P("## 5. One-at-a-time screen (1 run per cell)", "")
B = {s: BM.loc[s, "mean"] for s in ORDER}; NB = {s: int(BM.loc[s, "size"]) for s in ORDER}
c0r = oat[oat.oat_case == 0].set_index("strategy").load
P("Noise check, OAT case 0 (single run) over the pooled baseline mean: " +
  "; ".join(f"{s} {c0r[s] / B[s]:.2f}" for s in ORDER) + ".", "")
def ranking(o):
    v = o.sort_values(); items = list(v.items()); txt = items[0][0]
    for (a, x), (b, y) in zip(items, items[1:]):
        tie = abs(np.log(y / x)) < 2 * np.sqrt(S[a] ** 2 + S[b] ** 2)
        txt += (" ~ " if tie else " < ") + b
    return txt
def canon(txt):   # order-insensitive within ties
    return " < ".join(" ~ ".join(sorted(p.split(" ~ "))) for p in txt.split(" < "))
base_rank = canon(ranking(pd.Series(B)))
P(f"Tie rule for two single runs: |ln ratio| < 2 x sqrt(s_a^2 + s_b^2) with the pooled SDs of section 2 "
  f"(biological vs hybrid: within a factor {np.exp(2 * np.sqrt(S['biological'] ** 2 + S['hybrid'] ** 2)):.2f}). "
  f"Baseline ranking (pooled means): {base_rank}.", "")
rows, R = [], {}
for c in range(1, 19):
    o = oat[oat.oat_case == c].set_index("strategy").load
    rk = canon(ranking(o)); R[c] = {s: o[s] / B[s] for s in ORDER}
    flag = {s: abs(np.log(R[c][s])) > 2 * S[s] * np.sqrt(1 + 1 / NB[s]) for s in ORDER}
    rows.append([c, OAT[c]] + [f"{R[c][s]:.2f}{'' if flag[s] else ' (ns)'}" for s in ORDER] + [rk, "no" if rk == base_rank else "YES"])
P("Load ratio to the pooled baseline; (ns) = within 2 SD of replicate noise.", "")
P(*tbl(pd.DataFrame(rows, columns=["case", "assumption"] + ORDER + ["ranking", "changes?"])))
n_same = sum(r[-1] == "no" for r in rows)
P(f"Ranking unchanged in {n_same} of 18 cases. 'none' is last in 18 of 18. "
  f"Chemical (n_half 20) is worse than biological in {sum(oat[oat.oat_case == c].set_index('strategy').load['chemical'] > oat[oat.oat_case == c].set_index('strategy').load['biological'] for c in range(1, 19))} of 18.", "")
inf = pd.DataFrame({p: {s: max(abs(np.log10(R[lo][s])), abs(np.log10(R[hi][s]))) for s in ORDER} for p, (lo, hi) in PARAMS.items()}).T
inf["max_control"] = inf[["chemical", "biological", "hybrid"]].max(axis=1)
inf = inf.sort_values("max_control", ascending=False)
P("Influence = largest |log10 load ratio| over the low and high case (1.0 = tenfold).", "")
P(*tbl(inf.reset_index().rename(columns={"index": "assumption"}), {c: "{:.2f}".format for c in ORDER + ["max_control"]}))

# ---------- 6 Secondary outputs
P("## 6. Secondary outputs (main comparison means)", "")
none_load = t[t.strategy == "none"]["mean"].iloc[0]
s6 = t[["config", "ever", "end", "chem", "bio", "mean"]].copy()
s6["apps_per_farm_yr"] = (s6.chem + s6.bio) / N_FARMS / YEARS
s6["averted_per_app"] = (none_load - s6["mean"]) / (s6.chem + s6.bio).replace(0, np.nan)
s6["reduction_vs_none"] = 1 - s6["mean"] / none_load
P(f"Farms = {N_FARMS:,} farmer agents; applications are farm-level events; per farm per year averages over ALL farms, treated or not.", "")
P(*tbl(s6[["config", "ever", "end", "apps_per_farm_yr", "averted_per_app", "reduction_vs_none"]],
       {"ever": "{:.3f}".format, "end": "{:.3f}".format, "apps_per_farm_yr": "{:.2f}".format, "averted_per_app": "{:,.0f}".format, "reduction_vs_none": "{:.2%}".format}))
cy = main.groupby("strategy").cane.agg(["min", "max"])
P(f"Cane per harvested ha: none {cy.loc['none', 'min']:.2f}-{cy.loc['none', 'max']:.2f}; all controlled configurations "
  f"{main[main.strategy != 'none'].cane.min():.2f}-{main[main.strategy != 'none'].cane.max():.2f} (Y0 = 61.64; set by the assumed damage slope).", "")

# ---------- 7 Checks
P("## 7. Checks", "")
rows = []
for s in ORDER:
    a, b = wind[(wind.oat_case == 135) & (wind.strategy == s)].load, wind[(wind.oat_case == 315) & (wind.strategy == s)].load
    r = b.mean() / a.mean()
    rows.append([s, a.mean(), b.mean(), r, float(np.exp(2 * S[s])), "no difference" if abs(np.log(r)) < 2 * S[s] else "differs"])
P("Wind direction at bias 1.0 (n = 2 each): ratio 315/135 vs noise band exp(2 x s).", "")
P(*tbl(pd.DataFrame(rows, columns=["strategy", "load 135", "load 315", "ratio", "band", "verdict"]),
       {"load 135": "{:,.0f}".format, "load 315": "{:,.0f}".format, "ratio": "{:.2f}".format, "band": "{:.2f}".format}))
sc = sat.groupby("strategy").agg(load=("load", "mean"), cane=("cane", "mean")).loc[ORDER]
lc = b_main.groupby("strategy").agg(load=("load", "mean"), cane=("cane", "mean")).loc[ORDER]
P("Saturating vs linear damage (n = 2 each). Pest load cannot differ by construction (damage does not act on the pest); only cane can.", "")
P(*tbl(pd.DataFrame({"strategy": ORDER, "load linear": lc.load.values, "load saturating": sc.load.values, "cane linear": lc.cane.values,
                     "cane saturating": sc.cane.values}),
       {"load linear": "{:,.0f}".format, "load saturating": "{:,.0f}".format, "cane linear": "{:.2f}".format, "cane saturating": "{:.2f}".format}))

# ---------- 8 Scale contrast with the earlier island runs (random start; not pooled)
isl = os.path.join(A, "..", "main", "batch_main.csv")
if os.path.exists(isl):
    I = pd.read_csv(isl, header=None, names=COLS); I = I[I.experiment == "main"]
    im = I.groupby(["strategy", "n_half"]).load.agg(["size", "mean"])
    P("## 8. Earlier island runs (random start over 1% of cells, superseded design; context only)", "")
    P(*tbl(im.reset_index(), {"mean": "{:,.0f}".format}))
    try:
        P(f"none / chemical (no decline): island random start {im.loc[('none', 20), 'mean'] / im.loc[('chemical', 0), 'mean']:.0f}x; "
          f"window A seeded {none_load / c0:.0f}x. Different start and extent: not comparable beyond the sign.", "")
    except KeyError:
        pass

# ---------- 9 Extras (15 Nov): decline break-even, hybrid test, island hindcast probes
dec_p, hyb_p = os.path.join(A, "main", "batch_decline.csv"), os.path.join(A, "main", "batch_hybridtest.csv")
dec = None
if os.path.exists(dec_p):
    dec = pd.read_csv(dec_p, header=None, names=COLS)
    dm = dec.groupby("n_half").load.agg(["size", "mean", "min", "max"]).sort_index(ascending=False)
    dm["ratio_to_bio_hybrid"] = dm["mean"] / target
    P("## 9. Chemical at slower decline rates (window A, 2 runs each)", "",
      f"Biological/hybrid pool: mean {target:,.0f}, range {bh.min():,.0f}-{bh.max():,.0f}.", "")
    P(*tbl(dm.reset_index(), {"mean": "{:,.0f}".format, "min": "{:,.0f}".format, "max": "{:,.0f}".format, "ratio_to_bio_hybrid": "{:.2f}".format}))
if os.path.exists(hyb_p):
    ht = pd.read_csv(hyb_p, header=None, names=COLS)
    hm2 = ht.groupby(["oat_case", "strategy"]).agg(n=("load", "size"), mean=("load", "mean"), lo=("load", "min"), hi=("load", "max"),
                                                   chem=("apps_chem", "mean"), bio=("apps_bio", "mean")).reset_index()
    P("## 10. Hybrid test: chemical WITHOUT decline, conditions where farms reach Level 3 (window A, 2 runs each)", "",
      "Case 14 = treat only at Level 3; case 16 = re-treat every 28 days.", "")
    P(*tbl(hm2, {c: "{:,.0f}".format for c in ["mean", "lo", "hi", "chem", "bio"]}))
    for c in (14, 16):
        q = hm2[hm2.oat_case == c].set_index("strategy")["mean"]
        P(f"Case {c}: hybrid/biological {q['hybrid'] / q['biological']:.2f}; hybrid/chemical {q['hybrid'] / q['chemical']:.2f}; "
          f"biological/chemical {q['biological'] / q['chemical']:.2f}. Noise band for 2-run means: about "
          f"{np.exp(2 * np.sqrt(S['biological'] ** 2 / 2 + S['hybrid'] ** 2 / 2)):.2f}.", "")
SD = os.path.join(A, "..", "seeded")
HC = [("hybrid, baseline", "series", "hybrid_baseline"), ("none, baseline", "series", "none_baseline"),
      ("none, long jump 0.2", "hindcast", "none_j02"), ("none, transmission x4", "hindcast", "none_b4"),
      ("none, both", "hindcast", "none_b4j02"), ("none, transmission x2.5", "hindcast", "none_b25"),
      ("none, transmission x3", "hindcast", "none_b30")]
N_ISL, SRA_VAL, SRA_REP = 356869, 72103.42, 122009.81   # mapped cane cells; FOI Table I validated and reported ha (Occ + Or), 2 Oct 2026
SEED = ["Silay City", "City of Talisay", "City of Victorias", "Cadiz City", "Enrique B. Magalona (Saravia)", "Manapla", "Toboso",
        "Bago City", "La Carlota City", "La Castellana", "Murcia"]
hrows = []
for name, d, tag in HC:
    sp, lp = os.path.join(SD, d, f"series_{tag}.csv"), os.path.join(SD, d, f"lgu_{tag}_d731.csv")
    if not (os.path.exists(sp) and os.path.exists(lp)): continue
    L = pd.read_csv(lp, header=None, names=["tag", "seed", "idx", "name", "cane", "ever", "valid"])
    a = L.groupby("name").agg(cane=("cane", "first"), ever=("ever", "mean"), valid=("valid", "first")); a = a[(a.cane > 0)]
    ha = a.ever.sum(); v = a[a.valid >= 0]; ns_ = v[~v.index.isin(SEED)]
    rho = (v.ever / v.cane).rank().corr((v.valid / v.cane).rank()); rho2 = (ns_.ever / ns_.cane).rank().corr((ns_.valid / ns_.cane).rank())
    hrows.append([name, L.seed.count() // L.name.nunique(), ha, ha / SRA_VAL, int((a.ever > 0).sum()), rho, rho2])
if hrows:
    HT = pd.DataFrame(hrows, columns=["run", "n", "ha ever infested, 2 Oct 2026", "ratio to SRA validated", "LGUs infested", "Spearman", "Spearman non-seed"])
    P("## 11. Island hindcast, model day 731 = 2 Oct 2026 (seeded 22 May 2025)", "",
      f"SRA FOI Table I, 2 Oct 2026: validated {SRA_VAL:,.0f} ha, reported {SRA_REP:,.0f} ha. 1 cell = 1 ha; {N_ISL:,} mapped cane cells.", "")
    P(*tbl(HT, {"ha ever infested, 2 Oct 2026": "{:,.0f}".format, "ratio to SRA validated": "{:.2f}".format, "Spearman": "{:.2f}".format, "Spearman non-seed": "{:.2f}".format}))
    hd = HT.set_index("run")["ha ever infested, 2 Oct 2026"]
    if {"none, transmission x3", "none, transmission x4"} <= set(hd.index):
        f = np.log(SRA_VAL / hd["none, transmission x3"]) / np.log(hd["none, transmission x4"] / hd["none, transmission x3"])
        P(f"Tuning (no control): transmission x3 gives {hd['none, transmission x3'] / SRA_VAL:.2f} of the validated area; log-log interpolation between x3 and x4 "
          f"puts the match at about x{3 * (4 / 3) ** f:.2f} [inference; single runs]. One free parameter fitted to one number, with control ignored: tuned, not calibrated.", "")

# ---------- 12 Trajectory: model weekly curves vs dated official area figures
import datetime as _dt
OBS = [  # (date, ha, metric, scope, source)
    ("2025-05-22", 87.04, "affected, at declaration (= the model seed)", "Negros", "SRA via Inquirer / Watchmen 2025"),
    ("2025-06-19", 1505, "infested (SRA statement; date approximate)", "Negros Island Region", "Watchmen 20 Jun 2025"),
    ("2025-07-24", 2876.28, "infested (Res. 771)", "Negros Occidental only", "Manila Bulletin 25 Jul 2025"),
    ("2025-09-11", 4923, "affected (SRA status report)", "Negros Island", "USDA RP2025-0043 Table 2"),
    ("2026-06-26", 16019, "validated", "Negros Island", "PNA 29 Jun 2026 (saved capture)"),
    ("2026-08-14", 58374, "validated", "Negros Island", "SRA table via Digicast 1 Sep 2026 (secondary)"),
    ("2026-10-02", 72103.42, "validated (tuning / hindcast point)", "Negros Island", "SRA FOI Table I"),
]
D0 = _dt.date(2025, 5, 22)
def mday(d): return 233 + (_dt.date.fromisoformat(d) - D0).days
CUR = {}
for name, d, tag in HC:
    sp = os.path.join(SD, d, f"series_{tag}.csv")
    if os.path.exists(sp):
        q = pd.read_csv(sp, header=None); CUR[name] = q.groupby(5)[7].mean() * N_ISL
TUN = os.path.join(SD, "hindcast", "series_none_tuned.csv")
if os.path.exists(TUN):
    q = pd.read_csv(TUN, header=None); CUR["none, tuned transmission"] = q.groupby(5)[7].mean() * N_ISL
if CUR:
    rows = []
    for d, ha, metric, scope, src in OBS:
        md_ = mday(d); wk = md_ - md_ % 7
        rows.append([d, md_, f"{ha:,.0f}", metric + "; " + scope] + ["87 (seed)" if md_ == 233 else (f"{CUR[k].loc[wk]:,.0f}" if wk in CUR[k].index else "-") for k in CUR])
    P("## 12. Trajectory: model area ever infested (ha, island) vs dated official figures", "",
      "Model value = weekly series at the week containing the date. Observed metrics and scopes differ between rows (see column); "
      "early figures reflect detection catching up and SRA attributes part of the mid-2026 rise to validation effort. Sources are in OBS in analysis.py; the manuscript cites only the 22 May 2025, 11 Sep 2025 and 2 Oct 2026 figures.", "")
    P(*tbl(pd.DataFrame(rows, columns=["date", "model day", "observed ha", "metric; scope"] + list(CUR))))
    def growth(a, b, ser=None):
        return None
    o = {d: ha for d, ha, *_ in OBS}
    P("Observed multiplication between dates vs model (no control):", "")
    pr = [("2025-05-22", "2025-09-11"), ("2025-09-11", "2026-06-26"), ("2026-06-26", "2026-10-02")]
    rr = []
    for a, b in pr:
        row = [f"{a} to {b}", (mday(b) - mday(a)), f"{o[b] / o[a]:.1f}x"]
        for k in CUR:
            wa, wb = mday(a) - mday(a) % 7, mday(b) - mday(b) % 7
            c = CUR[k]
            if a == "2025-05-22": va = 87.0
            else: va = c.loc[wa] if wa in c.index else np.nan
            vb = c.loc[wb] if wb in c.index else np.nan
            row.append(f"{vb / va:.1f}x" if va and not np.isnan(vb) else "-")
        rr.append(row)
    P(*tbl(pd.DataFrame(rr, columns=["period", "days", "observed"] + list(CUR))))

# ---------- 13 Symmetry test: biological agent also loses efficacy (window A, 2 runs each; hybrid chemical arm at n_half 20)
bd_p = os.path.join(A, "main", "batch_biodecline.csv")
if os.path.exists(bd_p):
    bd = pd.read_csv(bd_p, header=None, names=COLS)
    bd["n_half_bio"] = bd.experiment.str.replace("biodec", "").astype(int)
    q = bd.groupby(["n_half_bio", "strategy"]).agg(n=("load", "size"), mean=("load", "mean"), lo=("load", "min"), hi=("load", "max"),
                                                   chem=("apps_chem", "mean"), bio=("apps_bio", "mean")).reset_index()
    chem_ref = {20: c20}
    if dec is not None and 80 in set(dec.n_half): chem_ref[80] = dec[dec.n_half == 80].load.mean()
    q["vs chemical, same decline"] = [m / chem_ref[n] if n in chem_ref else np.nan for n, m in zip(q.n_half_bio, q["mean"])]
    q["vs biological, no decline"] = q["mean"] / bm.mean()
    P("## 13. Symmetry test: the biological agent also halves in efficacy every n_half_bio applications per farm", "",
      f"References: chemical n_half 20 = {c20:,.0f}; chemical n_half 80 = {chem_ref.get(80, float('nan')):,.0f}; biological without decline (n = {len(bm)}) = {bm.mean():,.0f}. "
      "Hybrid keeps separate application counts for its two agents.", "")
    P(*tbl(q, {c: "{:,.0f}".format for c in ["mean", "lo", "hi", "chem", "bio"]} | {"vs chemical, same decline": "{:.2f}".format, "vs biological, no decline": "{:.1f}".format}))
    for n in sorted(set(q.n_half_bio)):
        w = q[q.n_half_bio == n].set_index("strategy")["mean"]
        P(f"n_half_bio {n}: hybrid/biological {w['hybrid'] / w['biological']:.2f}" + (f"; hybrid/chemical {w['hybrid'] / chem_ref[n]:.2f}" if n in chem_ref else "") + ".", "")

open(os.path.join(HERE, "numbers.md"), "w", encoding="utf-8").write("\n".join(out) + "\n")

# ---------- Figures
plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})
# (a) main comparison with replicate points
fig, ax = plt.subplots(figsize=(9, 4.8))
for i, r in t.iterrows():
    d = main[(main.strategy == r.strategy) & (main.n_half == r.n_half)]
    ax.plot(d.load, [i] * len(d), "o", ms=9, color=CB[r.strategy], alpha=0.85, mec="white")
    ax.plot([r["mean"]], [i], "|", ms=22, mew=2, color=CB[r.strategy])
    ax.text(r.hi * 1.25, i, f"{r['mean']:,.0f}", va="center", fontsize=10, color="#444")
ax.set_yticks(range(len(t))); ax.set_yticklabels(t.config); ax.invert_yaxis(); ax.set_xscale("log")
ax.set_xlim(t.lo.min() / 1.6, t.hi.max() * 6)
ax.set_xlabel("Cumulative pest load over 1,095 days (K-days, log scale)")
ax.set_title("Window A, seeded start: each dot is one run, bar = mean of 2", fontsize=12, loc="left")
ax.grid(axis="x", color="#ddd", lw=0.6); plt.tight_layout(); fig.savefig(os.path.join(HERE, "fig-a-main.png"), dpi=200); plt.close(fig)

# (b) OAT tornado
fig, axs = plt.subplots(1, 4, figsize=(13, 5), sharey=True, sharex=True)
names = list(inf.index)[::-1]
for ax, s in zip(axs, ORDER):
    for y, p in enumerate(names):
        lo, hi = PARAMS[p]
        ax.barh(y + 0.19, np.log10(R[hi][s]), height=0.36, color=CB[s])
        ax.barh(y - 0.19, np.log10(R[lo][s]), height=0.36, color=CB[s], alpha=0.45)
    band = np.log10(np.exp(2 * S[s] * np.sqrt(1 + 1 / NB[s])))
    ax.axvspan(-band, band, color="#999", alpha=0.25, lw=0); ax.axvline(0, color="k", lw=0.8)
    ax.set_title(s + (" (n_half 20)" if s in ("chemical", "hybrid") else ""), fontsize=12)
    ax.set_xlabel("log10 load ratio to baseline"); ax.grid(axis="x", color="#ddd", lw=0.6)
axs[0].set_yticks(range(len(names))); axs[0].set_yticklabels(names)
fig.suptitle("One-at-a-time screen, 1 run per bar: dark = high value, pale = low value; grey band = replicate noise (2 SD)", fontsize=11, x=0.01, ha="left")
plt.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(HERE, "fig-b-oat-tornado.png"), dpi=200); plt.close(fig)

# (c) chemical load vs efficacy decline
fig, ax = plt.subplots(figsize=(8, 4.8))
ch = main[main.strategy == "chemical"][["n_half", "load"]]
if dec is not None: ch = pd.concat([ch, dec[["n_half", "load"]]])
order = [n for n in [0, 160, 80, 40, 20, 10] if n in set(ch.n_half)]
pos = {n: i for i, n in enumerate(order)}
ax.axhspan(bh.min(), bh.max(), color=CB["biological"], alpha=0.25, lw=0)
ax.text(len(order) - 0.6, target, f"biological and hybrid\n(n = {len(bh)} runs)", va="center", ha="right", fontsize=10, color="#00644a")
ax.plot(ch.n_half.map(pos), ch.load, "o", ms=9, color=CB["chemical"], mec="white", zorder=3)
m = ch.groupby("n_half").load.mean(); ax.plot([pos[n] for n in order], [m[n] for n in order], "-", color=CB["chemical"], lw=1.5)
ax.axhline(none_load, color="k", ls="--", lw=1); ax.text(0.0, none_load * 1.15, "no control", fontsize=10)
ax.set_yscale("log"); ax.set_xticks(range(len(order))); ax.set_xticklabels(["no\ndecline" if n == 0 else str(n) for n in order])
ax.set_ylabel("Cumulative pest load (K-days, log scale)"); ax.set_xlabel("Chemical efficacy halves every ... applications per farm")
ax.set_title("Chemical control vs efficacy decline, window A (each dot = 1 run)", fontsize=12, loc="left")
ax.set_ylim(bh.min() / 2.2, none_load * 2.2); ax.grid(axis="y", color="#ddd", lw=0.6)
plt.tight_layout(); fig.savefig(os.path.join(HERE, "fig-c-decline.png"), dpi=200); plt.close(fig)

# (d) island hindcast
if hrows:
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    cols = [CB["hybrid"]] + ["#555"] * (len(HT) - 1)
    ax.barh(range(len(HT)), HT.iloc[:, 2], color=cols)
    for i, vv in enumerate(HT.iloc[:, 2]): ax.text(vv + 3000, i, f"{vv:,.0f}", va="center", fontsize=10, bbox=dict(fc="white", ec="none", pad=1.5), zorder=5)
    ax.axvline(SRA_VAL, color=CB["chemical"], lw=2); ax.text(SRA_VAL + 2000, -0.5, "SRA validated\n72,103 ha", color=CB["chemical"], fontsize=10, va="top")
    ax.axvline(SRA_REP, color=CB["chemical"], lw=1.5, ls="--"); ax.text(SRA_REP + 2000, -0.5, "SRA reported\n122,010 ha", color=CB["chemical"], fontsize=10, va="top")
    ax.set_yticks(range(len(HT))); ax.set_yticklabels(HT.run); ax.invert_yaxis(); ax.set_ylim(len(HT) - 0.2, -0.6)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlabel("Model cane area ever infested on 2 Oct 2026 (ha), all of Negros")
    ax.set_title("Hindcast from the 22 May 2025 seed: model vs SRA record", fontsize=12, loc="left")
    plt.tight_layout(); fig.savefig(os.path.join(HERE, "fig-d-hindcast.png"), dpi=200); plt.close(fig)
# (e) trajectory
if CUR:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sty = {"hybrid, baseline": (CB["hybrid"], "-"), "none, baseline": ("#000", "-"), "none, long jump 0.2": ("#777", "--"),
           "none, transmission x4": ("#777", "-"), "none, both": ("#bbb", "-"), "none, transmission x3": (CB["biological"], "-"), "none, transmission x2.5": (CB["biological"], ":")}
    for k, c in CUR.items():
        c = c[c.index >= 238]
        st = sty.get(k, (CB["biological"], "-"))
        ax.plot(c.index - 233, c.values, color=st[0], ls=st[1], lw=2, label="model: " + k + (" (tuned to 2 Oct 2026)" if k.endswith("x3") else ""))
    ax.plot([mday(d) - 233 for d, *_ in OBS], [h for _, h, *_ in OBS], "o", ms=9, color=CB["chemical"], mec="white", zorder=5, label="official figures (mixed metrics)")
    ax.set_yscale("log"); ax.set_xlim(0, 520); ax.set_ylim(50, 6e5)
    ax.set_xlabel("Days since 22 May 2025 (SRA declaration; model seed)"); ax.set_ylabel("Cane area ever infested (ha, log scale)")
    ax.set_title("Outbreak trajectory on Negros: model curves vs dated official figures", fontsize=12, loc="left")
    ax.legend(fontsize=8.5, frameon=False, loc="lower right"); ax.grid(axis="y", color="#ddd", lw=0.6)
    plt.tight_layout(); fig.savefig(os.path.join(HERE, "fig-e-trajectory.png"), dpi=200); plt.close(fig)
print("wrote numbers.md and figures in", HERE)
