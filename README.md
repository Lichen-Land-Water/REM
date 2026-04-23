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

## Repository purpose

This workflow is intended for terrain-based stream and floodplain analysis where you want to:

- derive stream centerlines from a DEM
- place transects perpendicular to a centerline
- identify low-elevation points that approximate a channel or water-surface control
- interpolate a continuous water surface
- compute a relative elevation raster
- convert the REM to mapped classes for interpretation or design

## Workflow overview

### Step 0: Extract streams and add attributes

Script: `0_get_streams.py`

Main function:

```python
def get_streams(
    dem: str,
    output_dir: str,
    threshold_km2: float = 1.0,
    overwrite: bool = False,
    breach_depressions: bool = True,
    thin_n: int = 10,
    create_thinned: bool = True,
    precip_raster: Optional[str] = None,
)
```

What it does:

- validates that the DEM is in a projected CRS
- breaches or fills depressions
- computes D8 flow direction and flow accumulation
- converts a drainage area threshold in km² to contributing cells
- extracts streams to raster and vector
- optionally creates a thinned centerline layer
- adds the following attributes to the stream GeoPackage:
  - `DA_km2`
  - `BF_width_Legg_m`
  - `BF_depth_Legg_m`
  - `BF_width_Castro_m`
  - `BF_depth_Castro_m`
  - `BF_width_Beechie_m`
  - `BF_depth_Beechie_m`

Notes:

- The drainage-area threshold is converted from km² to cell count using raster resolution and projected units.
- The script warns or fails if the DEM is not projected.
- The current Legg precipitation term is a placeholder basin-average value, not a raster-derived precipitation surface.

Typical output:

- stream raster
- stream shapefile
- stream GeoPackage
- optional thinned GeoPackage
- intermediate hydrologic rasters such as breached DEM, D8 pointer, and flow accumulation

### Step 1A: Create fixed-length transects

Script: `1_0_create_transects.py`

Main function:

```python
def create_bendy_transects_smooth(
    input_gpkg: str,
    output_gpkg: str,
    input_layer: Optional[str] = None,
    output_layer: Optional[str] = None,
    spacing: float = 100.0,
    transect_length: float = 1000.0,
    window: float = 200.0,
) -> str:
```

What it does:

- reads centerlines from a GeoPackage
- repairs invalid geometries where possible
- projects internally to a local UTM CRS if the input is geographic
- sorts centerlines by descending `DA_km2`
- places transects at a specified spacing
- computes a smoothed local perpendicular direction using a forward/backward window
- attempts to avoid overlapping transects by bending the line when needed
- writes transects back to the source CRS

Key output attributes include:

- `station`
- `centerline_id`
- `DA_km2`
- bankfull width/depth fields copied from the centerline layer

Use this version when you want one constant transect length across the network.

### Step 1B: Create drainage-area-scaled transects

Script: `1_1_create_transects_DA_weighted.py`

Main function:

```python
def create_bendy_transects_smooth(
    input_gpkg: str,
    output_gpkg: str,
    input_layer: Optional[str] = None,
    output_layer: Optional[str] = None,
    spacing: float = 100.0,
    window: float = 200.0,
) -> str:
```

What it does:

This version is similar to `1_0_create_transects.py`, but transect length varies by drainage area. The implemented formula is:

```python
transect_length = (DA_km2 ** (1/3)) * 200.0
```

It also writes an additional field:

- `transect_length_m`

Use this version when you want transects to scale with channel size.

### Step 2: Sample elevations along transects

Script: `2_0_get_elevations_along_transect.py`

Main functions:

```python
def extract_elevations_along_transect(
    transect_gpkg: str,
    dem_path: str,
    output_gpkg: str,
    flank_min_points: bool = False,
    method: str = "min",
)
```

```python
def extract_min_points(
    transect_gpkg: str,
    dem_path: str,
    output_gpkg: str,
    flank_min_points: bool = False,
)
```

```python
def extract_median_points(
    transect_gpkg: str,
    dem_path: str,
    output_gpkg: str,
    layer_name: str = "median_elev_points",
    flank_points: bool = True,
)
```

What it does:

- samples the DEM along each transect at approximately raster-resolution spacing
- supports `LineString` and `MultiLineString`
- identifies either:
  - the minimum-elevation sample point, or
  - the sampled point closest to the transect median elevation
- optionally also writes the start and end vertices of each transect
- carries through drainage area and bankfull fields where present

Important behavior:

- `method="min"` extracts minimum-elevation points
- `method="median"` extracts median-elevation control points
- `flank_min_points=True` includes the transect endpoints in addition to the central point
- null and non-positive elevations are dropped before writing output in the minimum-point workflow

Typical output:

- point GeoPackage with an `elevation` field

### Step 3: Interpolate water surface and compute REM

Script: `3_interpolate_water_surface.py`

Main functions:

```python
def interpolate_water_surface(
    gpkg_path: str,
    out_path: str,
    field: str,
    pix_size: float,
    power: float,
    smoothing: float,
    radius: float = None
) -> None:
```

```python
def difference_rasters(
    raster_path1: str,
    raster_path2: str,
    output_path: str,
    resampling: Resampling = Resampling.bilinear,
    out_dtype=np.float32,
    nodata_out=None,
)
```

```python
def merge_tifs(input_folder: str, output_path: str) -> None:
```

What it does:

- interpolates a raster water surface from point elevations using GDAL Grid with inverse distance weighting
- filters out points where the chosen interpolation field is null
- allows control over pixel size, IDW power, smoothing, and search radius
- computes a difference raster on the DEM grid:

```text
REM = DEM - interpolated_water_surface
```

Notes:

- the interpolation field can be `elevation` or another numeric field, such as a bankfull-depth attribute
- the difference raster is written on the grid of `raster_path1`
- `merge_tifs` is available if you need to mosaic tiled outputs

Typical output:

- interpolated water-surface raster
- REM raster

### Step 4A: Classify REM by bankfull-depth proportion

Script: `4_0_classify_rem_by_BFD.py`

Main function:

```python
def classify_rem_by_bankfull(
    rem_raster_path: str,
    output_class_raster_path: str,
    output_polygons_path: str,
    *,
    bf_raster_path: Optional[str] = None,
    bf_static_value: Optional[float] = None,
    thresholds: Sequence[float] = (0.5, 1.0, 2.0),
    out_nodata: int = 0,
    polygon_driver: str = "GPKG",
    polygon_layer: str = "rem_bf_classes",
    dissolve_polygons: bool = True,
) -> None:
```

What it does:

- classifies REM values by ratio to bankfull depth
- accepts either:
  - a bankfull-depth raster, or
  - a single static bankfull-depth value
- computes:

```text
ratio = REM / BF
```

- assigns classes by threshold breaks
- leaves values greater than the last threshold as unclassified nodata
- polygonizes the classified raster to a GeoPackage

Example interpretation for default thresholds `(0.5, 1.0, 2.0)`:

- class 1: `ratio <= 0.5`
- class 2: `0.5 < ratio <= 1.0`
- class 3: `1.0 < ratio <= 2.0`
- `ratio > 2.0` becomes nodata/unclassified

Polygon output fields include:

- `class_id`
- `Proportion of BF stage`

### Step 4B: Classify REM by manual bins

Script: `4_1_classify_rem_manual.py`

Main function:

```python
def classify_raster_manual_bins(
    in_raster_path: str,
    output_class_raster_path: str,
    output_polygons_path: str,
    *,
    class_edges: Sequence[float],
    out_nodata: int = 0,
    polygon_driver: str = "GPKG",
    polygon_layer: str = "classes",
    dissolve_polygons: bool = True,
) -> None:
```

What it does:

- classifies a raster using manual bin edges
- creates one class between each pair of consecutive edges
- uses inclusive upper bound only for the last class
- leaves values outside the specified range as nodata/unclassified
- polygonizes the result

Polygon output fields include:

- `class_id`
- `ClassRange`

Use this version when you want explicit REM classes rather than bankfull-scaled classes.

## Suggested processing sequence

### Option A: Full workflow from DEM

1. Run `0_get_streams.py` to derive streams and drainage-area/bankfull fields.
2. Run either `1_0_create_transects.py` or `1_1_create_transects_DA_weighted.py`.
3. Run `2_0_get_elevations_along_transect.py` to generate control points.
4. Run `3_interpolate_water_surface.py` to build the interpolated surface and REM.
5. Run either `4_0_classify_rem_by_BFD.py` or `4_1_classify_rem_manual.py`.

### Option B: Start from an existing centerline

If you already have a centerline GeoPackage, you can skip Step 0 and begin with transect generation.

## Inputs and outputs

### Core inputs

- projected DEM in meters or feet
- stream or centerline GeoPackage for transect generation
- optional bankfull-depth raster or static bankfull-depth value for classification

### Core outputs

- stream centerline GeoPackage
- transect GeoPackage
- minimum- or median-elevation point GeoPackage
- interpolated water-surface raster
- REM raster
- classified REM raster
- classified polygon GeoPackage

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

A typical install may look like this:

```bash
pip install geopandas numpy shapely pyproj rasterio fiona rasterstats whitebox whitebox-workflows scikit-learn
```

GDAL installation is environment-specific. In many setups it is easiest to install through conda-forge.

## Example usage

### Extract streams

```python
from 0_get_streams import get_streams

get_streams(
    dem=r"C:\path\to\dem.tif",
    output_dir=r"C:\path\to\outputs",
    threshold_km2=3,
    overwrite=False,
    breach_depressions=True,
    create_thinned=False,
)
```

### Create fixed-length transects

```python
from 1_0_create_transects import create_bendy_transects_smooth

create_bendy_transects_smooth(
    input_gpkg=r"C:\path\to\centerline.gpkg",
    output_gpkg=r"C:\path\to\transects.gpkg",
    output_layer="transects",
    spacing=100,
    transect_length=300,
    window=500,
)
```

### Create drainage-area-scaled transects

```python
from 1_1_create_transects_DA_weighted import create_bendy_transects_smooth

create_bendy_transects_smooth(
    input_gpkg=r"C:\path\to\centerline.gpkg",
    output_gpkg=r"C:\path\to\transects.gpkg",
    output_layer="transects",
    spacing=300,
    window=1000,
)
```

### Extract minimum-elevation points

```python
from 2_0_get_elevations_along_transect import extract_elevations_along_transect

extract_elevations_along_transect(
    transect_gpkg=r"C:\path\to\transects.gpkg",
    dem_path=r"C:\path\to\dem.tif",
    output_gpkg=r"C:\path\to\min_points.gpkg",
    flank_min_points=False,
    method="min",
)
```

### Build interpolated water surface and REM

```python
from 3_interpolate_water_surface import interpolate_water_surface, difference_rasters

interpolate_water_surface(
    gpkg_path=r"C:\path\to\min_points.gpkg",
    out_path=r"C:\path\to\water_surface.tif",
    field="elevation",
    pix_size=3,
    power=2,
    smoothing=1,
    radius=300,
)

difference_rasters(
    raster_path1=r"C:\path\to\dem.tif",
    raster_path2=r"C:\path\to\water_surface.tif",
    output_path=r"C:\path\to\rem.tif",
)
```

### Classify by bankfull-depth proportion

```python
from 4_0_classify_rem_by_BFD import classify_rem_by_bankfull

classify_rem_by_bankfull(
    rem_raster_path=r"C:\path\to\rem.tif",
    output_class_raster_path=r"C:\path\to\rem_bf_classes.tif",
    output_polygons_path=r"C:\path\to\rem_bf_classes.gpkg",
    bf_static_value=2.22,
    thresholds=(0.5, 1.0, 2.0),
    polygon_layer="floodplain",
    dissolve_polygons=True,
)
```

### Classify by manual bins

```python
from 4_1_classify_rem_manual import classify_raster_manual_bins

classify_raster_manual_bins(
    in_raster_path=r"C:\path\to\rem.tif",
    output_class_raster_path=r"C:\path\to\rem_classes.tif",
    output_polygons_path=r"C:\path\to\rem_classes.gpkg",
    class_edges=[-10, 0, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    polygon_layer="floodplain",
    dissolve_polygons=False,
)
```

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
