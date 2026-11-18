"""Join FOI Table I (LGU names) to HDX COD-AB adm3 P-codes; attach GEDI-S2 cane area. Name matches are explicit, checked by hand."""
import sys, csv, re
base = sys.argv[1]; foi = sys.argv[2]
cane = {(r['province'], r['adm3_name']): r for r in csv.DictReader(open(f"{base}/derived/cane_by_lgu.csv"))}
def norm(s):
    s = s.upper().replace('CITY OF ', '').replace(' CITY', '').replace('(CAPITAL)', '')
    s = re.sub(r'\(.*?\)', '', s).replace('STA.', 'SANTA').replace('-', '').replace('.', '').replace(' ', '')
    return s
manual = {'DONSALVADORBENEDICTO': 'SALVADORBENEDICTO', 'EBMAGALONA': 'ENRIQUEBMAGALONA', 'HINOBAAN': 'HINOBAAN', 'GUIHULNGAN': 'GUIHULNGAN', 'ZAMBOANGITA': 'ZAMBOANGUITA'}
idx = {(p, norm(n)): (p, n) for (p, n) in cane}
out = []; miss = []
for r in csv.DictReader(open(foi)):
    if r['subgroup_by_position_under_file_subtotal'].upper() in ('TOTAL',) or r['subgroup_by_position_under_file_subtotal'].startswith('SUBTOTAL'): continue
    k = norm(r['lgu_or_particulars']); k = manual.get(k, k)
    hit = idx.get((r['province'], k))
    if not hit: miss.append(r['lgu_or_particulars']); continue
    c = cane[hit]
    out.append([r['province'], r['lgu_or_particulars'], c['adm3_name'], c['adm3_pcode'], r['validated_total_area_affected_ha'], r['total_treated_area_ha'], round(float(c['gedis2_cane_ha']), 1), round(float(c['lgu_area_ha']), 1)])
with open(f"{base}/derived/foi_t1_crosswalk.csv", 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['province', 'foi_lgu', 'adm3_name', 'adm3_pcode', 'validated_ha_2oct2026', 'treated_ha_2oct2026', 'gedis2_cane_ha', 'lgu_area_ha']); w.writerows(out)
print('matched', len(out), 'missing', miss)
used = {(o[0], o[2]) for o in out}
print('adm3 without FOI row:', [n for (p, n) in cane if (p, n) not in used])
for o in sorted(out, key=lambda x: -x[6]): print(o[3], o[2][:24].ljust(24), 'cane', o[6], 'valid', o[4], 'v/c', round(float(o[4]) / o[6], 2) if o[6] else '-')
