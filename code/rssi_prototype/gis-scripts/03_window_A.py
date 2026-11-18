"""Window A (Cadiz City, Manapla, City of Victorias; 16-ha islet dropped): 1-ha grids in UTM 51N for GAMA.
cane_frac.asc  share of 10 m cane pixels per 1-ha cell (0-1), -1 outside the three LGUs
lgu.asc        LGU index (1 Cadiz, 2 Manapla, 3 Victorias), 0 outside
Cane cell rule (cane_frac >= 0.5) is an assumption, applied in GAMA, not here."""
import sys, numpy as np, geopandas as gpd, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
base = sys.argv[1]; out = f"{base}/derived/window_A"
names = ['Cadiz City', 'Manapla', 'City of Victorias']
a = gpd.read_file(f"{base}/derived/negros_adm3_utm51n.shp"); a = a[a.adm3_name.isin(names)]
parts = a.explode(index_parts=False); parts = parts[parts.area > 1e6]  # drop the 16-ha islet
d = parts.dissolve()
minx, miny, maxx, maxy = d.total_bounds
x0, y1 = np.floor(minx / 100) * 100, np.ceil(maxy / 100) * 100
W, H = int(np.ceil((maxx - x0) / 100)), int(np.ceil((y1 - miny) / 100))
T = from_origin(x0, y1, 100, 100)
frac = np.zeros((H, W), 'float32')
with rasterio.open(f"{base}/derived/negros_gedis2_cane.tif") as s:
    src = (s.read(1) == 1).astype('float32')
    reproject(src, frac, src_transform=s.transform, src_crs=s.crs, dst_transform=T, dst_crs='EPSG:32651', resampling=Resampling.average)
lgu = np.zeros((H, W), 'int16')
for i, n in enumerate(names, 1):
    lgu = np.where(rasterize(parts[parts.adm3_name == n].geometry, (H, W), transform=T, all_touched=False) == 1, i, lgu)
frac = np.where(lgu > 0, frac, -1).round(3)
def asc(p, arr, nd, fmt):
    with open(p, 'w') as f:
        f.write(f"ncols {W}\nnrows {H}\nxllcorner {x0}\nyllcorner {y1 - H*100}\ncellsize 100\nNODATA_value {nd}\n")
        np.savetxt(f, arr, fmt=fmt)
asc(f"{out}/cane_frac.asc", frac, -1, '%.3f'); asc(f"{out}/lgu.asc", lgu, 0, '%d')
for p in ('cane_frac.prj', 'lgu.prj'):
    open(f"{out}/{p}", 'w').write(rasterio.crs.CRS.from_epsg(32651).to_wkt())
parts.to_file(f"{out}/window_A_lgus.shp")
inside = lgu > 0; cane = (frac >= 0.5) & inside
print('grid', W, 'x', H, '=', W*H, '| inside LGUs', int(inside.sum()), '| cane cells (frac>=0.5)', int(cane.sum()), '| sum frac ha', round(float(frac[inside].sum()), 1))
for i, n in enumerate(names, 1): print(n, 'cells', int((lgu == i).sum()), 'cane cells', int(((lgu == i) & cane).sum()), 'frac ha', round(float(frac[lgu == i].sum()), 1))
