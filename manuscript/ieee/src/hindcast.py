import sys,csv,collections
from scipy.stats import spearmanr
R=sys.argv[1]
SEEDS={"Silay City","City of Talisay","City of Victorias","Cadiz City","Enrique B. Magalona (Saravia)","Manapla","Toboso","Bago City","La Carlota City","La Castellana","Murcia"}
for f in ["series/lgu_hybrid_baseline_d731","series/lgu_none_baseline_d731","hindcast/lgu_none_j02_d731","hindcast/lgu_none_b25_d731","hindcast/lgu_none_b30_d731","hindcast/lgu_none_b4_d731","hindcast/lgu_none_b4j02_d731"]:
    rows=[r for r in csv.reader(open(f"{R}/seeded/{f}.csv")) if r]
    runs=collections.defaultdict(list)
    for r in rows: runs[r[1]].append(r)
    tot=[sum(int(r[5]) for r in v) for v in runs.values()]
    # per LGU mean over runs
    d=collections.defaultdict(lambda:[0,0,0,None])
    for v in runs.values():
        for r in v:
            d[r[3]][0]=int(r[4]); d[r[3]][1]+=int(r[5])/len(runs); d[r[3]][2]=float(r[6])
    a=[(k,e/c,v/c) for k,(c,e,v,_) in d.items() if c>0 and v>=0]
    rho=spearmanr([x[1] for x in a],[x[2] for x in a])[0]
    b=[x for x in a if x[0] not in SEEDS]
    rho2=spearmanr([x[1] for x in b],[x[2] for x in b])[0]
    print(f, len(runs), round(sum(tot)/len(tot)), round(sum(tot)/len(tot)/72103,2), len(a), round(rho,2), len(b), round(rho2,2))
