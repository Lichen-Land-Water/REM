# get_streams.py

import math
import os
import warnings
from typing import Optional

import fiona
import numpy as np
import geopandas as gpd
import rasterio
import whitebox
from rasterstats import zonal_stats
from shapely.geometry import LineString, MultiLineString
from whitebox_workflows import WbEnvironment


M2_TO_SQMI = 3.861021585424458e-7
SQMI_TO_KM2 = 2.589988110336
IN_TO_CM = 2.54
MM_TO_IN = 1.0 / 25.4
FT_TO_M = 0.3048


def _get_raster_area_metadata(raster_path: str) -> dict:
    """
    Returns raster CRS and pixel area information needed to convert
    cell counts to physical area.

    Returns
    -------
    dict with keys:
        crs
        pixel_width
        pixel_height
        pixel_area_native
        pixel_area_m2
        linear_units
    """
    with rasterio.open(raster_path) as src:
        crs = src.crs
        xres, yres = src.res
        pixel_width = abs(xres)
        pixel_height = abs(yres)

    if crs is None:
        raise ValueError(f"Raster has no CRS: {raster_path}")

    if not crs.is_projected:
        raise ValueError(
            f"Raster is not projected: {raster_path} ({crs}). "
            "Drainage area and threshold conversions based on pixel area are not reliable "
            "for geographic CRS. Reproject the DEM/raster to a projected CRS first."
        )

    try:
        linear_units = crs.linear_units.lower() if crs.linear_units else None
    except Exception:
        linear_units = None

    pixel_area_native = pixel_width * pixel_height

    if linear_units in {"metre", "meter", "metres", "meters", "m"}:
        pixel_area_m2 = pixel_area_native
    elif linear_units in {"foot", "feet", "us survey foot", "foot_us", "ft"}:
        pixel_area_m2 = pixel_area_native * 0.09290304
    else:
        raise ValueError(
            f"Unrecognized projected CRS linear units '{linear_units}' for raster: {raster_path}"
        )

    return {
        "crs": crs,
        "pixel_width": pixel_width,
        "pixel_height": pixel_height,
        "pixel_area_native": pixel_area_native,
        "pixel_area_m2": pixel_area_m2,
        "linear_units": linear_units,
    }


def km2_to_cell_threshold(flow_accum_raster: str, threshold_km2: float) -> int:
    """
    Convert a drainage-area threshold in km² to a flow-accumulation threshold
    in number of contributing cells.
    """
    meta = _get_raster_area_metadata(flow_accum_raster)
    threshold_m2 = threshold_km2 * 1_000_000.0
    threshold_cells = math.ceil(threshold_m2 / meta["pixel_area_m2"])
    return int(threshold_cells)


def add_raster_mean_to_streams(
    streams_gpkg_path: str,
    raster_path: str,
    output_field: str,
    all_touched: bool = True,
) -> None:
    """
    Adds the mean raster value intersecting each stream feature to a GeoPackage.

    CRS handling
    ------------
    - The stream GeoPackage keeps its original CRS.
    - A temporary copy of the stream geometry is reprojected to the raster CRS
      before raster values are assigned.
    """
    if raster_path is None:
        warnings.warn(f"No raster path provided for {output_field}. Skipping.")
        return

    streams = gpd.read_file(streams_gpkg_path)
    streams = streams[streams.geometry.notnull()].copy()
    streams = streams[streams.is_valid].copy()

    if streams.empty:
        warnings.warn(
            f"No valid stream features found in {streams_gpkg_path}. "
            f"Skipping {output_field}."
        )
        return

    if streams.crs is None:
        raise ValueError(
            f"Streams layer has no CRS in {streams_gpkg_path}. "
            "Set the stream CRS before assigning raster values."
        )

    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
        nodata = src.nodata

    if raster_crs is None:
        raise ValueError(f"Raster has no CRS: {raster_path}")

    if streams.crs != raster_crs:
        print(
            f"Reprojecting temporary stream geometry from {streams.crs} "
            f"to raster CRS {raster_crs} for {output_field} extraction..."
        )
        streams_for_stats = streams.to_crs(raster_crs)
    else:
        streams_for_stats = streams.copy()

    stats = zonal_stats(
        vectors=streams_for_stats.geometry,
        raster=raster_path,
        stats=["mean"],
        nodata=nodata,
        all_touched=all_touched,
    )

    streams[output_field] = [
        s["mean"] if s.get("mean") is not None else None
        for s in stats
    ]

    streams.to_file(streams_gpkg_path, driver="GPKG")


def add_PRISM_to_streams(
    streams_gpkg_path: str,
    ppt_raster_path: str,
    tmean_raster_path: str,
) -> None:
    """
    Adds PRISM 30-year average precipitation and mean temperature to streams.

    Stored output fields
    --------------------
    - ann_precip_in: annual precipitation in inches
    - tmean_degC: mean annual temperature in degrees C

    Expected input raster units
    ---------------------------
    - ppt raster: mm/year
    - tmean raster: degrees C
    """
    add_raster_mean_to_streams(
        streams_gpkg_path=streams_gpkg_path,
        raster_path=ppt_raster_path,
        output_field="__ppt_mm_yr",
    )

    add_raster_mean_to_streams(
        streams_gpkg_path=streams_gpkg_path,
        raster_path=tmean_raster_path,
        output_field="tmean_degC",
    )

    streams = gpd.read_file(streams_gpkg_path)

    if "__ppt_mm_yr" in streams.columns:
        streams["ann_precip_in"] = streams["__ppt_mm_yr"] * MM_TO_IN
        streams = streams.drop(columns="__ppt_mm_yr")

    streams.to_file(streams_gpkg_path, driver="GPKG")


def get_streams(
    dem: str,
    output_dir: str,
    threshold_km2: float = 1.0,
    overwrite: bool = False,
    breach_depressions: bool = True,
    thin_n: int = 10,
    create_thinned: bool = True,
    precip_raster: Optional[str] = None,
    temp_raster: Optional[str] = None,
):
    """
    Process a DEM to extract streams, save them to a GeoPackage, optionally
    create a thinned centerline layer, and add drainage area, PRISM climate
    values, and bankfull dimensions.

    Stored stream-network fields are limited to:
    - DA_sqmi
    - ann_precip_in
    - tmean_degC

    Bankfull functions convert these values internally where equations require
    km² or cm.

    Parameters
    ----------
    dem : str
        Input DEM path.
    output_dir : str
        Output directory.
    threshold_km2 : float, default 1.0
        Stream initiation threshold expressed as drainage area in square kilometers.
    overwrite : bool, default False
        Whether to overwrite intermediate outputs.
    breach_depressions : bool, default True
        If True, breach depressions. If False, fill depressions.
    thin_n : int, default 10
        Keep every nth vertex if create_thinned is True.
    create_thinned : bool, default True
        Whether to create a thinned centerline layer.
    precip_raster : Optional[str], default None
        PRISM precipitation raster path, expected in mm/year.
    temp_raster : Optional[str], default None
        PRISM mean temperature raster path, expected in degrees C.

    Returns
    -------
    str
        Path to stream GeoPackage.
    """
    wbt = whitebox.WhiteboxTools()
    _ = WbEnvironment()

    try:
        dem_meta = _get_raster_area_metadata(dem)
    except ValueError as e:
        warnings.warn(str(e))
        raise

    os.makedirs(output_dir, exist_ok=True)

    filled_dem = os.path.join(output_dir, "filled_dem.tif")
    d8_pointer = os.path.join(output_dir, "d8_pointer.tif")
    flow_accum = os.path.join(output_dir, "flow_accum.tif")
    breached_dem = os.path.join(output_dir, "breached_dem.tif")

    if overwrite:
        for path in [filled_dem, d8_pointer, flow_accum, breached_dem]:
            if os.path.exists(path):
                os.remove(path)

    if breach_depressions and not os.path.exists(breached_dem):
        wbt.breach_depressions_least_cost(dem, breached_dem, 10)

    if not breach_depressions and not os.path.exists(filled_dem):
        wbt.fill_depressions(dem, filled_dem)

    src_dem = breached_dem if breach_depressions else filled_dem

    if not os.path.exists(d8_pointer) or not os.path.exists(flow_accum):
        wbt.d8_pointer(src_dem, d8_pointer)
        wbt.d8_flow_accumulation(src_dem, flow_accum)

    try:
        threshold_cells = km2_to_cell_threshold(flow_accum, threshold_km2)
    except ValueError as e:
        warnings.warn(str(e))
        raise

    threshold_label = str(threshold_km2).replace(".", "p")
    streams_raster = os.path.join(output_dir, f"streams_{threshold_label}km2.tif")

    if not os.path.exists(streams_raster):
        print(
            f"Extracting streams using threshold = {threshold_km2} km² "
            f"({threshold_cells} contributing cells)"
        )
        wbt.extract_streams(flow_accum, streams_raster, threshold_cells)

    streams_shp = streams_raster.replace(".tif", ".shp")
    streams_gpkg = streams_raster.replace(".tif", ".gpkg")
    streams_layer = os.path.splitext(os.path.basename(streams_gpkg))[0]

    if not os.path.exists(streams_shp):
        wbt.raster_streams_to_vector(streams_raster, d8_pointer, streams_shp)

    gdf = gpd.read_file(streams_shp)
    dem_crs = dem_meta["crs"]

    if gdf.crs is None:
        gdf = gdf.set_crs(dem_crs)
    else:
        gdf = gdf.to_crs(dem_crs)

    gdf.to_file(streams_gpkg, layer=streams_layer, driver="GPKG")

    thinned_gpkg = None
    if create_thinned:
        print(f"Thinning centerline by keeping every {thin_n}th vertex...")
        thinned_gpkg = streams_gpkg.replace(".gpkg", "_thinned.gpkg")
        thin_centerline(
            input_gpkg=streams_gpkg,
            layer_name=streams_layer,
            output_gpkg=thinned_gpkg,
            output_layer=f"{streams_layer}_thinned",
            n=thin_n,
        )

    print("Adding drainage area in square miles...")
    add_DA_to_stream(streams_gpkg, flow_accum)
    print("Adding DEM-derived longitudinal slope...")
    add_dem_slope_to_streams(streams_gpkg, dem)
    
    if precip_raster is not None and temp_raster is not None:
        print("Adding PRISM precipitation in inches and temperature in degrees C...")
        add_PRISM_to_streams(
            streams_gpkg_path=streams_gpkg,
            ppt_raster_path=precip_raster,
            tmean_raster_path=temp_raster,
        )
    else:
        warnings.warn(
            "PRISM precipitation and/or temperature raster not provided. "
            "Bankfull equations that use precipitation will fall back to default values."
        )

    print("Adding bankfull dimensions...")
    add_BF_to_streams_Legg(streams_gpkg)
    add_BF_to_streams_Castro(streams_gpkg)
    add_BF_to_streams_Beechie(streams_gpkg)

    print(f"[✔] Streams extracted to: {streams_gpkg}")
    return streams_gpkg


def add_DA_to_stream(
    streams_gpkg: str,
    flow_accum_raster: str,
    da_field: str = "DA_sqmi",
    buffer_cells: float = 0.75,
) -> None:
    """
    Adds drainage area in square miles to stream features in `streams_gpkg`
    using a flow-accumulation raster.

    Assumes the flow accumulation raster stores upslope contributing area as
    number of cells.

    Notes
    -----
    - Uses the maximum flow-accumulation value within a small buffer around
      each stream feature.
    - Pixel area is calculated automatically from the raster resolution and CRS.
    - Warns and exits if the raster is not projected.
    """
    streams = gpd.read_file(streams_gpkg)

    if streams.empty:
        warnings.warn(
            f"No stream features found in {streams_gpkg}. "
            "Skipping drainage area calculation."
        )
        return

    try:
        meta = _get_raster_area_metadata(flow_accum_raster)
    except ValueError as e:
        warnings.warn(str(e))
        return

    raster_crs = meta["crs"]

    with rasterio.open(flow_accum_raster) as src:
        nodata = src.nodata

    if streams.crs is None:
        warnings.warn(
            f"Streams layer has no CRS in {streams_gpkg}. "
            f"Assuming raster CRS: {raster_crs}."
        )
        streams = streams.set_crs(raster_crs)
    elif streams.crs != raster_crs:
        streams = streams.to_crs(raster_crs)

    buffer_dist = max(meta["pixel_width"], meta["pixel_height"]) * buffer_cells
    buffers = streams.geometry.buffer(buffer_dist)

    stats = zonal_stats(
        buffers,
        flow_accum_raster,
        stats=["max"],
        nodata=nodata,
        all_touched=True,
    )

    streams[da_field] = [
        (s["max"] * meta["pixel_area_m2"] * M2_TO_SQMI)
        if s.get("max") is not None
        else None
        for s in stats
    ]

    streams.to_file(streams_gpkg, driver="GPKG")

def add_dem_slope_to_streams(
    streams_gpkg: str,
    dem_raster: str,
    slope_field: str = "slope_m_m",
    slope_pct_field: Optional[str] = "slope_pct",
    sample_spacing: Optional[float] = None,
    min_samples: int = 3,
    z_units_to_xy_units: float = 1.0,
) -> None:
    """
    Adds DEM-derived longitudinal slope to each stream feature.

    Method
    ------
    For each LineString or MultiLineString feature:
    1. Reproject stream geometry to the DEM CRS, if needed.
    2. Sample DEM elevations along the stream at regular streamwise spacing.
    3. Fit a least-squares line to elevation versus streamwise distance.
    4. Store the absolute value of the fitted slope.

    This is more appropriate for stream-network slope than assigning mean values
    from a raster slope grid, because it estimates slope along the flow path.

    Stored output fields
    --------------------
    - slope_m_m: dimensionless longitudinal slope
    - slope_pct: slope percent, optional

    Parameters
    ----------
    streams_gpkg : str
        Stream GeoPackage path.
    dem_raster : str
        DEM raster path.
    slope_field : str, default "slope_m_m"
        Output field for dimensionless slope.
    slope_pct_field : Optional[str], default "slope_pct"
        Optional output field for percent slope. Set to None to skip.
    sample_spacing : Optional[float], default None
        Spacing between DEM samples in DEM CRS units. If None, uses the smaller
        DEM cell dimension.
    min_samples : int, default 3
        Minimum valid DEM samples required to compute slope.
    z_units_to_xy_units : float, default 1.0
        Multiplier to convert DEM elevation units to DEM horizontal CRS units.
        Use 1.0 when horizontal and vertical units are the same.
        Example: if DEM elevations are meters but CRS units are feet, use 3.28084.
    """
    streams = gpd.read_file(streams_gpkg)

    if streams.empty:
        warnings.warn(
            f"No stream features found in {streams_gpkg}. "
            "Skipping DEM slope calculation."
        )
        return

    if streams.crs is None:
        raise ValueError(
            f"Streams layer has no CRS in {streams_gpkg}. "
            "Set the stream CRS before calculating DEM-derived slope."
        )

    with rasterio.open(dem_raster) as src:
        dem_crs = src.crs
        nodata = src.nodata
        xres, yres = src.res
        default_spacing = min(abs(xres), abs(yres))

        if dem_crs is None:
            raise ValueError(f"DEM has no CRS: {dem_raster}")

        spacing = sample_spacing or default_spacing

        if spacing <= 0:
            raise ValueError("sample_spacing must be positive.")

        if streams.crs != dem_crs:
            streams_for_slope = streams.to_crs(dem_crs)
        else:
            streams_for_slope = streams.copy()

        def _sample_line_profile(line: LineString, start_offset: float = 0.0):
            if line is None or line.is_empty or line.length <= 0:
                return [], []

            distances = np.arange(0.0, line.length, spacing).tolist()

            if not distances or distances[-1] < line.length:
                distances.append(line.length)

            points = [line.interpolate(d) for d in distances]
            coords = [(pt.x, pt.y) for pt in points]

            sampled = list(src.sample(coords))

            valid_distances = []
            valid_elevations = []

            for d, value_array in zip(distances, sampled):
                if value_array is None or len(value_array) == 0:
                    continue

                z = float(value_array[0])

                if nodata is not None and np.isclose(z, nodata):
                    continue

                if not np.isfinite(z):
                    continue

                valid_distances.append(start_offset + d)
                valid_elevations.append(z * z_units_to_xy_units)

            return valid_distances, valid_elevations

        def _profile_slope(geom):
            if geom is None or geom.is_empty:
                return None

            all_distances = []
            all_elevations = []
            offset = 0.0

            if isinstance(geom, LineString):
                distances, elevations = _sample_line_profile(geom, start_offset=0.0)
                all_distances.extend(distances)
                all_elevations.extend(elevations)

            elif isinstance(geom, MultiLineString):
                for part in geom.geoms:
                    distances, elevations = _sample_line_profile(part, start_offset=offset)
                    all_distances.extend(distances)
                    all_elevations.extend(elevations)
                    offset += part.length

            else:
                return None

            if len(all_distances) < min_samples:
                return None

            x = np.asarray(all_distances, dtype=float)
            z = np.asarray(all_elevations, dtype=float)

            valid = np.isfinite(x) & np.isfinite(z)

            if valid.sum() < min_samples:
                return None

            x = x[valid]
            z = z[valid]

            if np.nanmax(x) == np.nanmin(x):
                return None

            # Linear regression: z = a * distance + b.
            # Slope is reported as absolute longitudinal fall per unit length.
            a, _b = np.polyfit(x, z, 1)

            return abs(float(a))

        streams[slope_field] = streams_for_slope.geometry.apply(_profile_slope)

    if slope_pct_field is not None:
        streams[slope_pct_field] = streams[slope_field] * 100.0

    streams.to_file(streams_gpkg, driver="GPKG")

def _get_precip_in(streams: gpd.GeoDataFrame, fallback_precip_in: float = 72.17 / 2.54):
    """
    Returns annual precipitation in inches.

    The fallback is 72.17 cm converted to inches, consistent with the previous
    default precipitation value.
    """
    if "ann_precip_in" not in streams.columns:
        warnings.warn(
            "ann_precip_in not found. "
            f"Using fallback precipitation value of {fallback_precip_in:.2f} inches."
        )
        return fallback_precip_in

    return streams["ann_precip_in"]


def add_BF_to_streams_Legg(streams_gpkg_path: str) -> str:
    """
    Computes bankfull width and depth using Legg & Olson 2015.

    Stored input fields used
    ------------------------
    - DA_sqmi
    - ann_precip_in

    Internal conversions
    --------------------
    - DA_sqmi to DA_km2 for the depth equation
    - ann_precip_in to ann_precip_cm for the depth equation
    """
    streams = gpd.read_file(streams_gpkg_path)
    streams = streams[streams.geometry.notnull()].copy()
    streams = streams[streams.is_valid].copy()

    if "DA_sqmi" not in streams.columns:
        raise ValueError("DA_sqmi field is required for Legg bankfull equations.")

    da_sqmi = streams["DA_sqmi"]
    da_km2 = da_sqmi * SQMI_TO_KM2

    precip_in = _get_precip_in(streams)
    precip_cm = precip_in * IN_TO_CM

    streams["BF_width_Legg_m"] = (
        FT_TO_M
        * 1.16
        * 0.91
        * (da_sqmi ** 0.381)
        * (precip_in ** 0.634)
    )

    streams["BF_depth_Legg_m"] = (
        0.0939
        * (da_km2 ** 0.233)
        * (precip_cm ** 0.264)
    )

    streams.to_file(streams_gpkg_path, driver="GPKG")
    return streams_gpkg_path


def add_BF_to_streams_Castro(streams_gpkg_path: str) -> str:
    """
    Computes bankfull width and depth using Castro & Jackson 2001.

    Stored input fields used
    ------------------------
    - DA_sqmi

    Castro & Jackson equations use drainage area in square miles and return
    width/depth in feet, which are converted to meters.
    """
    streams = gpd.read_file(streams_gpkg_path)

    if "DA_sqmi" not in streams.columns:
        raise ValueError("DA_sqmi field is required for Castro bankfull equations.")

    da_sqmi = streams["DA_sqmi"]

    streams["BF_width_Castro_m"] = (
        FT_TO_M
        * 9.40
        * (da_sqmi ** 0.42)
    )

    streams["BF_depth_Castro_m"] = (
        FT_TO_M
        * 0.61
        * (da_sqmi ** 0.33)
    )

    streams.to_file(streams_gpkg_path, driver="GPKG")
    return streams_gpkg_path


def add_BF_to_streams_Beechie(streams_gpkg_path: str) -> str:
    """
    Computes bankfull width and depth using Beechie and Imaki 2013.

    Stored input fields used
    ------------------------
    - DA_sqmi
    - ann_precip_in

    Internal conversions
    --------------------
    - DA_sqmi to DA_km2
    - ann_precip_in to ann_precip_cm
    """
    streams = gpd.read_file(streams_gpkg_path)

    if "DA_sqmi" not in streams.columns:
        raise ValueError("DA_sqmi field is required for Beechie bankfull equations.")

    if "BF_width_Castro_m" not in streams.columns or "BF_depth_Castro_m" not in streams.columns:
        raise ValueError(
            "BF_width_Castro_m and BF_depth_Castro_m are required before computing "
            "Beechie depth. Run add_BF_to_streams_Castro first."
        )

    da_km2 = streams["DA_sqmi"] * SQMI_TO_KM2

    precip_in = _get_precip_in(streams)
    precip_cm = precip_in * IN_TO_CM

    streams["BF_width_Beechie_m"] = (
        0.177
        * (da_km2 ** 0.397)
        * (precip_cm ** 0.453)
    )

    streams["BF_depth_Beechie_m"] = (
        streams["BF_width_Beechie_m"]
        * streams["BF_depth_Castro_m"]
        / streams["BF_width_Castro_m"]
    )

    streams.to_file(streams_gpkg_path, driver="GPKG")
    return streams_gpkg_path


def thin_centerline(
    input_gpkg: str,
    layer_name: str,
    output_gpkg: str,
    output_layer: Optional[str] = None,
    n: int = 10,
) -> None:
    """
    Keeps every nth vertex in each LineString or part of a MultiLineString
    from the specified layer in input_gpkg, writing to output_gpkg.
    """
    gdf = gpd.read_file(input_gpkg, layer=layer_name)

    def _thin(geom):
        if geom is None or geom.is_empty:
            return geom

        if isinstance(geom, LineString):
            coords = list(geom.coords)
            pts = coords[::n]

            if coords[-1] not in pts:
                pts.append(coords[-1])

            return LineString(pts)

        if isinstance(geom, MultiLineString):
            parts = []

            for part in geom.geoms:
                coords = list(part.coords)
                pts = coords[::n]

                if coords[-1] not in pts:
                    pts.append(coords[-1])

                parts.append(LineString(pts))

            return MultiLineString(parts)

        return geom

    gdf["geometry"] = gdf.geometry.apply(_thin)

    out_layer = output_layer or layer_name
    gdf.to_file(output_gpkg, layer=out_layer, driver="GPKG")


def threshold_lines_by_length(
    input_gpkg: str,
    output_gpkg: str,
    threshold: float = 1200.0,
) -> None:
    """
    Read the first layer of `input_gpkg`, keep only LineString/MultiLineString
    features longer than `threshold` in CRS units, and write them to `output_gpkg`.
    """
    layers = fiona.listlayers(input_gpkg)

    if not layers:
        raise ValueError(f"No layers found in {input_gpkg!r}")

    layer = layers[0]
    gdf = gpd.read_file(input_gpkg, layer=layer)

    is_line = gdf.geometry.type.isin(["LineString", "MultiLineString"])
    lines = gdf[is_line].copy()

    lines["__length"] = lines.geometry.length
    filtered = lines[lines["__length"] > threshold].drop(columns="__length")

    filtered.to_file(output_gpkg, layer=layer, driver="GPKG")


if __name__ == "__main__":
    dem = r"C:\L\Lichen\Lichen - Documents\Projects\20260003_Owens-Snipe Assessment (UCSWCD)\07_GIS\0_Data_In\Public\LiDAR\USGS1m_proj_2020-2021_road_crossing_conditioned.tiff"
    threshold_km2 = 1
    output_dir = r"C:\L\Lichen\Lichen - Documents\Projects\20260003_Owens-Snipe Assessment (UCSWCD)\07_GIS\1_Analysis\Stream Network Analysis\Streams\1m"

    get_streams(
        dem=dem,
        output_dir=output_dir,
        threshold_km2=threshold_km2,
        overwrite=False,
        breach_depressions=True,
        create_thinned=False,
        precip_raster=r"C:\L\Lichen\Lichen - Documents\Library\GIS\PRISM\prism_ppt_30yr_avg_mmyr.tif",
        temp_raster=r"C:\L\Lichen\Lichen - Documents\Library\GIS\PRISM\prism_tmean_30yr_avg_degC.tif",
    )