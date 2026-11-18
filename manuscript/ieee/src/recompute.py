"""Recompute the result numbers in the IEEE paper from the raw batch CSVs. Usage: python recompute.py <results folder>"""
import sys, csv, statistics as st, math, collections
R=sys.argv[1]
C=["experiment","oat_case","strategy","n_half","seed","share_ever","share_end","load","chem","bio","cane","harv","maxD"]
def rd(p):
    return [dict(zip(C,r)) for r in csv.reader(open(f"{R}/{p}")) if r]
def g(rows,**k):
    return [r for r in rows if all(str(r[a])==str(b) for a,b in k.items())]
def m(rows,f="load"): return st.mean(float(r[f]) for r in rows)
out={}
A=rd("A/main/batch_main.csv")
for s,n in [("none","20"),("chemical","0"),("chemical","20"),("chemical","10"),("biological","20"),("hybrid","0"),("hybrid","20"),("hybrid","10")]:
    rows=g(A,strategy=s,n_half=n)
    if not rows: rows=g(A,strategy=s)
    out[f"A {s} {n}"]=(len(rows),round(m(rows)),round(m(rows,"share_ever"),3),round(m(rows,"chem")),round(m(rows,"bio")),round(sum(float(r['cane'])/float(r['harv']) for r in rows)/len(rows),2))
I=rd("seeded/main/batch_main.csv")
for s,n in [("none","20"),("chemical","0"),("chemical","20"),("chemical","10"),("biological","20"),("hybrid","20")]:
    rows=g(I,strategy=s,n_half=n)
    out[f"I {s} {n}"]=(len(rows),round(m(rows)),round(m(rows,"share_ever"),3),round(m(rows,"chem")),round(m(rows,"bio")))
D=rd("A/main/batch_decline.csv")
for n in ["40","80","160"]: out[f"decl {n}"]=(len(g(D,n_half=n)),round(m(g(D,n_half=n))))
H=rd("A/main/batch_hybridtest.csv")
for c in ["14","16"]:
    for s in ["chemical","biological","hybrid"]:
        r=g(H,oat_case=c,strategy=s); out[f"hyb {c} {s}"]=(len(r),round(m(r)))
B=rd("A/main/batch_biodecline.csv")
for r in B: pass
bd=collections.defaultdict(list)
for r in B: bd[(r['experiment'],r['strategy'])].append(float(r['load']))  # experiment = biodec20 / biodec80
for k,v in bd.items(): out[f"bio-decl {k}"]=(len(v),round(st.mean(v)))
S=rd("A/main/batch_saturating.csv"); O=rd("A/oat/batch_oat.csv")
pool={}
for s in ["none","chemical","biological","hybrid"]:
    v=[float(r['load']) for r in g(A,strategy=s,n_half="20")+g(S,strategy=s)+g(O,oat_case="0",strategy=s)]
    pool[s]=st.mean(v); out[f"pool {s}"]=(len(v),round(st.mean(v)))
bh=[float(r['load']) for r in g(A,strategy="biological")+g(S,strategy="biological")+g(O,oat_case="0",strategy="biological")+g(A,strategy="hybrid")+g(S,strategy="hybrid")+g(O,oat_case="0",strategy="hybrid")]
out["bio+hyb pool"]=(len(bh),round(st.mean(bh)),round(min(bh)),round(max(bh)))
P=st.mean(bh)
for k in ["A chemical 0","A chemical 20","A chemical 10"]: out["ratio "+k]=round(out[k][1]/P,2)
for n in ["40","80","160"]: out["ratio decl "+n]=round(out[f"decl {n}"][1]/P,2)
none=out["A none 20"][1]
for k in [k for k in out if k.startswith("A ") and "none" not in k]: out["red "+k]=round(100*(1-out[k][1]/none),2)
# OAT ranking changes
for c in range(1,19):
    loads={s:float(g(O,oat_case=str(c),strategy=s)[0]['load']) for s in ["none","chemical","biological","hybrid"]}
    out[f"oat {c}"]={s:round(loads[s]/pool[s],2) for s in loads}
for f in ["lgu_hybrid_baseline_d731","lgu_none_baseline_d731"]:
    pass
for k,v in out.items(): print(k,v)
