"""Island-wide 1-ha landscape for GAMA (all of Negros: 57 LGUs, both provinces).
Step 1: cane fraction in ~100 m blocks of 10x10 native pixels (EPSG:4326), read in strips.
Step 2: average-resample that fraction to 100 m UTM 51N cells.
Output derived/negros/: cane_frac_100m.tif, lgu_idx_100m.tif, lgu_index.csv"""
import sys, csv, numpy as np, geopandas as gpd, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin, Affine
base = sys.argv[1]; out = f"{base}/derived/negros"
src_p = f"{base}/derived/negros_gedis2_cane.tif"
with rasterio.open(src_p) as s:
    H, W = s.height, s.width; T = s.transform
    h2, w2 = H // 10, W // 10
    frac = np.zeros((h2, w2), 'float32')
    for i in range(h2):
        blk = s.read(1, window=((i * 10, i * 10 + 10), (0, w2 * 10))) == 1
        frac[i] = blk.reshape(10, w2, 10).mean(axis=(0, 2))
    T2 = T * Affine.scale(10, 10)
adm = gpd.read_file(f"{base}/derived/negros_adm3_utm51n.shp")
minx, miny, maxx, maxy = adm.total_bounds
x0, y1 = np.floor(minx / 100) * 100 - 200, np.ceil(maxy / 100) * 100 + 200
Wd, Hd = int(np.ceil((maxx + 200 - x0) / 100)), int(np.ceil((y1 - (miny - 200)) / 100))
TD = from_origin(x0, y1, 100, 100)
dst = np.zeros((Hd, Wd), 'float32')
reproject(frac, dst, src_transform=T2, src_crs='EPSG:4326', dst_transform=TD, dst_crs='EPSG:32651', resampling=Resampling.average)
adm = adm.sort_values('adm3_pcode').reset_index(drop=True); adm['idx'] = np.arange(1, len(adm) + 1)
lgu = rasterize(((g, i) for g, i in zip(adm.geometry, adm.idx)), (Hd, Wd), transform=TD, fill=0, dtype='int16')
prof = dict(driver='GTiff', width=Wd, height=Hd, count=1, crs='EPSG:32651', transform=TD, compress='deflate')
with rasterio.open(f"{out}/cane_frac_100m.tif", 'w', dtype='float32', **prof) as d: d.write(dst, 1)
with rasterio.open(f"{out}/lgu_idx_100m.tif", 'w', dtype='int16', **prof) as d: d.write(lgu, 1)
with open(f"{out}/lgu_index.csv", 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['idx', 'adm3_pcode', 'adm3_name', 'adm2_name'])
    for r in adm.itertuples(): w.writerow([r.idx, r.adm3_pcode, r.adm3_name, r.adm2_name])
inside = lgu > 0
print('grid', Wd, 'x', Hd, '=', Wd * Hd, '| inside LGUs', int(inside.sum()))
for th in (0.5,): print('cane cells (frac>=%.1f)' % th, int(((dst >= th) & inside).sum()), '| frac-sum ha inside', round(float(dst[inside].sum()), 1))
for p in ['Negros Occidental', 'Negros Oriental']:
    ids = adm.idx[adm.adm2_name == p].values; m = np.isin(lgu, ids)
    print(p, 'frac-sum ha', round(float(dst[m].sum()), 1), 'cane cells', int(((dst >= 0.5) & m).sum()))
