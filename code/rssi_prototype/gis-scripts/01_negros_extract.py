"""Extract Negros from HDX COD-AB (adm3) and the GEDI-S2 sugarcane map; zonal cane area per LGU.
Inputs (untouched): GIS/raw/hdx-codab/phl_admin_boundaries.gdb.zip, GIS/raw/zenodo-gedis2/philippines_GEDIS2_v1.zip
Outputs: GIS/derived/negros_adm3_utm51n.shp (EPSG:32651), negros_gedis2_cane.tif (band1 sugarcane 0/1, band2 tall months), cane_by_lgu.csv
"""
import sys, zipfile, numpy as np, geopandas as gpd, rasterio, pyogrio
from rasterio.merge import merge
from rasterio.mask import mask
base = sys.argv[1]
gdb = f"/vsizip/{base}/raw/hdx-codab/phl_admin_boundaries.gdb.zip/phl_admin_boundaries.gdb"
print(pyogrio.list_layers(gdb))
lyr = [l for l, _ in pyogrio.list_layers(gdb) if 'adm3' in l.lower() or 'admin3' in l.lower()][0]
adm3 = gpd.read_file(gdb, layer=lyr)
adm3 = adm3[adm3['adm2_name'].isin(['Negros Occidental', 'Negros Oriental'])].copy()
print(lyr, len(adm3), adm3.crs, adm3['valid_on'].unique(), adm3['version'].unique())
adm3[["adm3_name","adm3_pcode","adm2_name","adm2_pcode","geometry"]].to_crs(32651).to_file(f"{base}/derived/negros_adm3_utm51n.shp")
z = f"{base}/raw/zenodo-gedis2/philippines_GEDIS2_v1.zip"
tiles = [n for n in zipfile.ZipFile(z).namelist() if n.endswith('-0000058368.tif') and ('0000058368-' in n or '0000116736-' in n)]
xmin, ymin, xmax, ymax = adm3.to_crs(4326).total_bounds
srcs = [rasterio.open(f"/vsizip/{z}/{n}") for n in tiles]
arr, tr = merge(srcs, bounds=(xmin - .01, ymin - .01, xmax + .01, ymax + .01), indexes=[2, 1])
prof = srcs[0].profile.copy(); prof.update(height=arr.shape[1], width=arr.shape[2], transform=tr, count=2, compress='deflate', tiled=True, BIGTIFF='IF_SAFER')
out = f"{base}/derived/negros_gedis2_cane.tif"
with rasterio.open(out, 'w', **prof) as d:
    d.write(arr); d.set_band_description(1, 'sugarcane'); d.set_band_description(2, 'n_tallmonths')
print('raster', arr.shape, tr)
# pixel area by row latitude (m2), WGS84 approx
res = tr.a; R = 6371008.8
rows = np.arange(arr.shape[1]); lat = tr.f + (rows + .5) * tr.e
pix_ha = (np.radians(res) * R) ** 2 * np.cos(np.radians(lat)) / 1e4
recs = []
with rasterio.open(out) as d:
    for _, r in adm3.to_crs(4326).iterrows():
        m, mt = mask(d, [r.geometry], crop=True, indexes=1, filled=True, nodata=255)
        r0 = int(round((mt.f - tr.f) / tr.e))
        w = pix_ha[r0:r0 + m.shape[0]][:, None]
        recs.append((r['adm2_name'], r['adm3_name'], r['adm3_pcode'], float(((m == 1) * w).sum()), float(((m != 255) * w).sum())))
import csv
with open(f"{base}/derived/cane_by_lgu.csv", 'w', newline='') as f:
    c = csv.writer(f); c.writerow(['province', 'adm3_name', 'adm3_pcode', 'gedis2_cane_ha', 'lgu_area_ha']); c.writerows(recs)
for p in ['Negros Occidental', 'Negros Oriental']:
    print(p, round(sum(x[3] for x in recs if x[0] == p), 1), 'ha cane;', round(sum(x[4] for x in recs if x[0] == p), 1), 'ha land')
