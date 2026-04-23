# Relative Elevation Model (REM) Workflow

This repository contains a Python workflow for generating a stream-based Relative Elevation Model (REM) and optional classified floodplain products from a DEM and vectorized centerlines or streams.

The scripts are organized as a stepwise pipeline:

1. Extract streams from a DEM and add drainage area and bankfull attributes
2. Create cross-valley transects from stream centerlines
3. Sample elevations along those transects to identify low points
4. Interpolate a water surface from sampled points and subtract it from the DEM to create a REM
5. Classify the REM either by bankfull-depth proportions or by manual elevation bins

The repository currently includes the following scripts:

- `0_get_streams.py`
- `1_0_create_transects.py`
- `1_1_create_transects_DA_weighted.py`
- `2_0_get_elevations_along_transect.py`
- `3_interpolate_water_surface.py`
- `4_0_classify_rem_by_BFD.py`
- `4_1_classify_rem_manual.py`


## Dependencies

Based on the repository scripts, the workflow uses:

- `geopandas`
- `numpy`
- `shapely`
- `pyproj`
- `rasterio`
- `fiona`
- `rasterstats`
- `whitebox`
- `whitebox_workflows`
- `scikit-learn`
- `GDAL` / `osgeo`

Installation:

```bash
conda env create -f rem_environment.yml
```
Or if that doesn't work, try to let conda solve the environments itself:

```bash
conda install geopandas numpy shapely pyproj rasterio fiona rasterstats whitebox whitebox-workflows scikit-learn
```


## Workflow overview

### Step 0: Extract streams and add attributes

Script: `0_get_streams.py`

What it does:

- validates that the DEM is in a projected CRS
- breaches or fills depressions
- computes D8 flow direction and flow accumulation
- converts a drainage area threshold in km² to contributing cells
- extracts streams to raster and vector
- optionally creates a thinned centerline layer
- adds the following attributes to the stream GeoPackage:
  - `DA_km2`
  - 'ann_precip_mm'
  - `BF_width_Legg_m`
  - `BF_depth_Legg_m`
  - `BF_width_Castro_m`
  - `BF_depth_Castro_m`
  - `BF_width_Beechie_m`
  - `BF_depth_Beechie_m`

Notes:

- Streams may hae to be manually pruned after creation to select for reach of interest
- The drainage-area threshold is converted from km² to cell count using raster resolution and projected units.
- The script warns or fails if the DEM is not projected.

Outputs:
- stream GeoPackage
- intermediate hydrologic rasters such as breached DEM, D8 pointer, and flow accumulation

### Step 1: Create fixed-length transects

Script: `1_0_create_transects.py`

What it does:

- sorts centerlines by descending `DA_km2`
- places non intersecting transects at a specified spacing
- computes a smoothed local perpendicular direction using a forward/backward window
- attempts to avoid overlapping transects by bending the line when needed
- writes transects back to the source CRS

Output attributes include:

- `station`
- `centerline_id`
- `DA_km2`
- bankfull width/depth fields copied from the centerline layer

### Step 2: Sample elevations along transects

Script: `2_0_get_elevations_along_transect.py`

What it does:

- samples the DEM along each transect at specified spacing
- identifies either:
  - the minimum-elevation sample point (HAWS REM)
  - the sampled point closest to the transect median elevation (GGL REM)

Important behavior:

- `method="min"` extracts minimum-elevation points
- `method="median"` extracts median-elevation control points
- `flank_min_points=True` includes the transect endpoints in addition to the central point
- null and non-positive elevations are dropped before writing output in the minimum-point workflow

Typical output:

- point GeoPackage with an `elevation` field

### Step 3: Interpolate water surface and compute REM

Script: `3_interpolate_water_surface.py`

What it does:

- interpolates a raster water surface from point elevations using GDAL Grid with inverse distance weighting
- filters out points where the chosen interpolation field is null
- allows control over pixel size, IDW power, smoothing, and search radius
- computes a difference raster on the DEM grid:

```text
REM = DEM - interpolated_water_surface
```

Notes:

- the interpolation field can be `elevation` or another numeric field, such as a bankfull-depth attribute. This is useful for basin scale analysis of terrain height above bankfull depth

Typical output:

- interpolated water-surface raster
- REM raster


## Assumptions and caveats

- DEM and raster inputs should be in a projected CRS.
- Several operations assume linear units are meaningful for distance, area, and raster resolution.
- The drainage-area-scaled transect script currently uses `transect_length = (DA_km2 ** (1/3)) * 200.0`, even though one docstring line still mentions `DA_km2 * 5`. The implemented formula is the one actually used.
- The Legg precipitation term in `0_get_streams.py` is currently a fixed placeholder value.
- The interpolation approach is IDW-based and may require project-specific tuning of radius, smoothing, and pixel size.
- The minimum-elevation approach assumes the low point on each transect is a reasonable control for the interpolated surface. That may not hold in all geomorphic settings.

## File descriptions

### `0_get_streams.py`
Extract streams from a DEM, convert a drainage-area threshold from km² to contributing cells, vectorize streams, optionally thin centerlines, and append drainage-area and bankfull-geometry attributes.

### `1_0_create_transects.py`
Create fixed-length, smoothed, de-conflicted transects perpendicular to a centerline.

### `1_1_create_transects_DA_weighted.py`
Create drainage-area-scaled, smoothed, de-conflicted transects perpendicular to a centerline.

### `2_0_get_elevations_along_transect.py`
Sample the DEM along transects and export minimum- or median-elevation control points, with optional flank endpoints.

### `3_interpolate_water_surface.py`
Interpolate a water surface from point data using GDAL Grid and compute a REM by differencing rasters.

### `4_0_classify_rem_by_BFD.py`
Classify a REM by multiples of bankfull depth and polygonize the class raster.

### `4_1_classify_rem_manual.py`
Classify a REM using manual class edges and polygonize the class raster.
