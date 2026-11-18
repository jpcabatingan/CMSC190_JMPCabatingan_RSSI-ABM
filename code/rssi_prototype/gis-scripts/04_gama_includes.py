"""Model inputs for window A, written to sp-tracker/gama/rssi_prototype/includes/.
landscape_A.asc  value = lgu_id*10 + cane_frac (lgu 1 Cadiz, 2 Manapla, 3 Victorias); -9999 outside the 3 LGUs
lgus_A.shp       municipal boundaries (UTM 51N)
roads_A.shp      OSM trunk/primary/secondary/tertiary roads clipped to the grid extent (UTM 51N), field hwy
OSM data (c) OpenStreetMap contributors, ODbL."""
import sys, json, numpy as np, geopandas as gpd
from shapely.geometry import LineString, box
gis, inc = sys.argv[1], sys.argv[2]
w = f"{gis}/derived/window_A"
def read_asc(p):
    hdr = {}
    with open(p) as f:
        for _ in range(6):
            k, v = f.readline().split(); hdr[k.lower()] = float(v)
        a = np.loadtxt(f)
    return hdr, a
h, frac = read_asc(f"{w}/cane_frac.asc"); _, lgu = read_asc(f"{w}/lgu.asc")
val = np.where(lgu > 0, lgu * 10 + np.clip(frac, 0, 1), -9999)
with open(f"{inc}/landscape_A.asc", 'w') as f:
    f.write(f"ncols {int(h['ncols'])}\nnrows {int(h['nrows'])}\nxllcorner {h['xllcorner']}\nyllcorner {h['yllcorner']}\ncellsize 100\nNODATA_value -9999\n")
    np.savetxt(f, val, fmt='%.3f')
open(f"{inc}/landscape_A.prj", 'w').write(open(f"{w}/cane_frac.prj").read())
lg = gpd.read_file(f"{w}/window_A_lgus.shp")[['adm3_name', 'adm3_pcode', 'geometry']]
lg.to_file(f"{inc}/lgus_A.shp")
js = json.load(open(f"{gis}/raw/osm/roads_A_overpass.json"))
rows = [(e['tags'].get('highway'), e['tags'].get('name', ''), LineString([(p['lon'], p['lat']) for p in e['geometry']]))
        for e in js['elements'] if e['type'] == 'way' and len(e.get('geometry', [])) > 1]
r = gpd.GeoDataFrame({'hwy': [x[0] for x in rows], 'name': [x[1][:60] for x in rows]}, geometry=[x[2] for x in rows], crs=4326).to_crs(32651)
x0, y0 = h['xllcorner'], h['yllcorner']; x1, y1 = x0 + h['ncols'] * 100, y0 + h['nrows'] * 100
r = r.clip(box(x0, y0, x1, y1)); r = r[r.geometry.length > 0].explode(index_parts=False)
r = r[r.geom_type == 'LineString']
r.to_file(f"{inc}/roads_A.shp")
inside = val > -9999
print('cells', val.size, 'inside', int(inside.sum()), 'cane>=0.5', int(((val % 10) >= 0.5)[inside].sum()))
print(r.groupby('hwy').geometry.apply(lambda g: round(g.length.sum() / 1000, 1)).to_dict(), 'km')
