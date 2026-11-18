# GIS inputs for the Negros landscape (8 Nov 2026)

raw/ (downloaded 8 Nov 2026, unmodified; do not edit)
- zenodo-gedis2/philippines_GEDIS2_v1.zip: Di Tommaso, Wang, Strey & Lobell, sugarcane map 2019-2022, 10 m, EPSG:4326, v01 (25 Mar 2024), CC BY 4.0. https://zenodo.org/records/10871164 . Paper: ESSD 16:4931-4947 (2024), doi:10.5194/essd-16-4931-2024. Bands: n_tallmonths, sugarcane (0/1), ESA, ESRI, GLAD crop masks.
- hdx-codab/phl_admin_boundaries.gdb.zip, .xlsx: OCHA COD-AB Philippines, admin 0-4, valid_on 2025-02-13, version v03. Source PSA/NAMRIA. Licence CC BY-IGO. Indicative boundaries, not official. Negros Occidental still coded under Region VI (PH06045), Negros Oriental under Region VII (PH07046), not NIR.

scripts/
- 01_negros_extract.py: Negros adm3 to UTM 51N shapefile; Negros clip of the cane map; mapped cane ha per LGU.
- 02_foi_crosswalk.py: FOI Table I LGU names to adm3_pcode (56/56 matched; Bacong has no FOI row).

derived/
- negros_adm3_utm51n.shp (EPSG:32651), negros_gedis2_cane.tif (b1 sugarcane, b2 tall months), cane_by_lgu.csv, foi_t1_crosswalk.csv.

Check [computed]: mapped cane Negros Occidental 265,081.8 ha vs SRA CY2025-26 195,515.09 (x1.36); Negros Oriental 86,905.4 vs 43,700.53 (x1.99); Negros 351,987.2 vs 239,215.62 (x1.47). Same order of magnitude; map overstates. Validated area (2 Oct 2026) is below mapped cane in every LGU (max ratio 0.90, San Jose).
