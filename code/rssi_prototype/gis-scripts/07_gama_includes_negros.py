"""GAMA inputs for all of Negros, into sp-tracker/gama/rssi_prototype/includes/.
landscape_negros.asc  value = lgu_idx x 10 + cane fraction (lgu_idx 1-57, see lgu_index.csv); -9999 outside Negros LGUs
lgus_negros.shp       57 LGUs (UTM 51N): idx, adm3_name, adm3_pcode, adm2_name, valid_ha (SRA FOI Table I, 2 Oct 2026; -1 = no row)
roads_negros.shp      OSM trunk/primary/secondary inside the Negros LGUs (UTM 51N); (c) OpenStreetMap contributors, ODbL"""
import sys, csv, json, numpy as np, geopandas as gpd, rasterio
from shapely.geometry import LineString
gis, inc = sys.argv[1], sys.argv[2]; d = f"{gis}/derived/negros"
with rasterio.open(f"{d}/cane_frac_100m.tif") as s: frac = s.read(1); T = s.transform
with rasterio.open(f"{d}/lgu_idx_100m.tif") as s: lgu = s.read(1)
val = np.where(lgu > 0, lgu * 10 + np.clip(frac, 0, 1), -9999)
H, W = val.shape
with open(f"{inc}/landscape_negros.asc", 'w') as f:
    f.write(f"ncols {W}\nnrows {H}\nxllcorner {T.c}\nyllcorner {T.f - H * 100}\ncellsize 100\nNODATA_value -9999\n")
    np.savetxt(f, val, fmt='%.3f')
open(f"{inc}/landscape_negros.prj", 'w').write(rasterio.crs.CRS.from_epsg(32651).to_wkt())
idx = {r['adm3_pcode']: int(r['idx']) for r in csv.DictReader(open(f"{d}/lgu_index.csv"))}
foi = {r['adm3_pcode']: float(r['validated_ha_2oct2026']) for r in csv.DictReader(open(f"{gis}/derived/foi_t1_crosswalk.csv"))}
adm = gpd.read_file(f"{gis}/derived/negros_adm3_utm51n.shp")
adm['idx'] = adm.adm3_pcode.map(idx); adm['valid_ha'] = adm.adm3_pcode.map(foi).fillna(-1.0)
adm[['idx', 'adm3_name', 'adm3_pcode', 'adm2_name', 'valid_ha', 'geometry']].to_file(f"{inc}/lgus_negros.shp")
js = json.load(open(f"{gis}/raw/osm/roads_negros_overpass.json"))
rows = [(e['tags'].get('highway'), LineString([(p['lon'], p['lat']) for p in e['geometry']]))
        for e in js['elements'] if e['type'] == 'way' and len(e.get('geometry', [])) > 1]
r = gpd.GeoDataFrame({'hwy': [x[0] for x in rows]}, geometry=[x[1] for x in rows], crs=4326).to_crs(32651)
r = r.clip(adm.union_all()).explode(index_parts=False); r = r[(r.geom_type == 'LineString') & (r.length > 0)]
r.to_file(f"{inc}/roads_negros.shp")
print('asc', W, 'x', H, '| roads km', r.assign(km=r.length / 1000).groupby('hwy').km.sum().round(1).to_dict(), '| segments', len(r))
