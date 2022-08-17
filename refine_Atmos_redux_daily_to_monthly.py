#!/usr/bin/env python

# This script contains the refineDiags that reduce data in time
# from daily to monthly values (e.g. min/ax over the period)
# It can accept any file and will only compute refineDiags in fields
# are present.

import argparse
import os
import netCDF4 as nc
import xarray as xr


CMOR_MISSING_VALUE = 1.0e20
extra_time_variables = ["time_bnds", "average_T1", "average_T2", "average_DT"]
do_not_encode_vars = ["nv", "grid_xt", "grid_yt", "time"]
grid_vars = ["grid_xt", "grid_yt", "ap", "b", "ap_bnds", "b_bnds", "lev", "lev_bnds"]
pkgname = "xrefineDiags"
scriptname = "refine_Atmos_redux_daily_to_monthly.py"
pro = f"{pkgname}/{scriptname}"


def run():
    """run the refineDiags"""

    # --- parse command line arguemnts
    args = parse_args()
    verbose = args.verbose

    # --- open file
    if verbose:
        print(f"{pro}: Opening input file {args.infile}")

    ds = xr.open_dataset(args.infile, decode_cf=False)

    # --- create output dataset
    if os.path.exists(args.outfile):
        if verbose:
            print(f"{pro}: Opening existing file {args.outfile}")

        refined = xr.open_dataset(args.outfile, decode_cf=False)
        refined.load()  # required!
    else:
        if verbose:
            print(f"{pro}: Creating new dataset")

        refined = xr.Dataset()

    # -- create time variables if they don't exist in refined dataset
    if "time" not in list(refined.variables):
        refined = populate_time_vars(refined, ds)
        print(refined)

    # -- compute min/max from daily to monthly
    if "t_ref_max" in list(ds.variables):
        refined["t_ref_max"] = max_over_month(ds["t_ref_max"], ds, refined)
    if "t_ref_min" in list(ds.variables):
        refined["t_ref_min"] = min_over_month(ds["t_ref_min"], ds, refined)

    if "u_ref_max" in list(ds.variables):
        refined["u_ref_max"] = max_over_month(ds["u_ref_max"], ds, refined)
    if "u_ref_min" in list(ds.variables):
        refined["u_ref_min"] = min_over_month(ds["u_ref_min"], ds, refined)

    if "rh_ref_max" in list(ds.variables):
        refined["rh_ref_max"] = max_over_month(ds["rh_ref_max"], ds, refined)
    if "rh_ref_min" in list(ds.variables):
        refined["rh_ref_min"] = min_over_month(ds["rh_ref_min"], ds, refined)

    # --- write dataset to file
    new_vars_output = len(list(refined.variables)) > 0

    if verbose and new_vars_output:
        print(
            f"{pro}: writting variables {list(refined.variables)} into refined file {args.outfile} "
        )
    elif verbose and not new_vars_output:
        print(f"{pro}: no variables created, not writting refined file")

    if new_vars_output:
        write_dataset(refined, ds, args)

    return None


def max_over_month(da, ds, refined):

    realtime = xr.decode_cf(ds)["time"]
    da_max = da.groupby(realtime.dt.month).mean(dim="time").rename({"month": "time"})
    da_max["time"] = refined["time"]
    da_max.attrs = da.attrs.copy()
    da_max.attrs["cell_methods"] = da_max.attrs["cell_methods"].replace(
        "max", "max within days time: mean over days"
    )

    return da_max


def min_over_month(da, ds, refined):

    realtime = xr.decode_cf(ds)["time"]
    da_min = da.groupby(realtime.dt.month).mean(dim="time").rename({"month": "time"})
    da_min["time"] = refined["time"]
    da_min.attrs = da.attrs.copy()
    da_min.attrs["cell_methods"] = da_min.attrs["cell_methods"].replace(
        "min", "min within days time: mean over days"
    )

    return da_min


def populate_time_vars(refined, ds):

    realtime = xr.decode_cf(ds)["time"]

    time_monthly = (
        ds["time"].groupby(realtime.dt.month).mean(dim="time").rename({"month": "time"})
    )
    refined["time"] = xr.DataArray(
        data=time_monthly, dims=("time"), attrs=ds["time"].attrs
    )

    average_DT = (
        ds["average_DT"]
        .groupby(realtime.dt.month)
        .sum(dim="time")
        .rename({"month": "time"})
    )
    refined["average_DT"] = xr.DataArray(data=average_DT.values, dims=("time"))
    refined["average_DT"].attrs = ds["average_DT"].attrs.copy()

    average_T1 = (
        ds["average_T1"]
        .groupby(realtime.dt.month)
        .min(dim="time")
        .rename({"month": "time"})
    )
    refined["average_T1"] = xr.DataArray(data=average_T1.values, dims=("time"))
    refined["average_T1"].attrs = ds["average_T1"].attrs.copy()

    average_T2 = (
        ds["average_T2"]
        .groupby(realtime.dt.month)
        .max(dim="time")
        .rename({"month": "time"})
    )
    refined["average_T2"] = xr.DataArray(data=average_T2.values, dims=("time"))
    refined["average_T2"].attrs = ds["average_T2"].attrs.copy()

    time_bnds = xr.concat([average_T1, average_T2], dim="nv").transpose()
    refined["time_bnds"] = xr.DataArray(data=time_bnds.values, dims=("time", "nv"))
    refined["time_bnds"].attrs = ds["time_bnds"].attrs.copy()

    # needs to be repeated to conserve attributes
    refined["time"] = xr.DataArray(
        data=time_monthly, dims=("time"), attrs=ds["time"].attrs
    )

    refined["nv"] = xr.DataArray(
        [1.0, 2.0], dims="nv", attrs={"long_name": "vertex number"}
    )

    return refined


def write_dataset(ds, template, args):
    """prepare the dataset and dump into netcdf file"""

    if len(ds.attrs) == 0:
        ds.attrs = template.attrs.copy()  # copy global attributes
    ds.attrs["filename"] = args.outfile

    # --- add proper grid attrs
    for var in grid_vars:
        if var in list(template.variables):
            ds[var] = template[var]
            ds[var].attrs = template[var].attrs.copy()

    # --- remove bounds in attributes since it messed the bnds var
    var_with_bounds = []
    bounds_variables = []
    for var in list(ds.variables):
        if "bounds" in ds[var].attrs:
            var_with_bounds.append(var)
            bounds_variables.append(ds[var].attrs.pop("bounds"))

    encoding = set_netcdf_encoding(ds)

    ds.to_netcdf(
        args.outfile, format=args.format, encoding=encoding, unlimited_dims="time"
    )
    post_write(args.outfile, var_with_bounds, bounds_variables)

    return None


def set_netcdf_encoding(ds):
    """set preferred options for netcdf encoding"""

    all_vars = list(ds.variables)
    encoding = {}

    for var in do_not_encode_vars:
        if var in all_vars:
            encoding.update({var: dict(_FillValue=None)})

    return encoding


def post_write(filename, var_with_bounds, bounds_variables):
    """fix a posteriori attributes that xarray.to_netcdf
    did not do properly using low level netcdf lib"""

    f = nc.Dataset(filename, "a")

    for var, bndvar in zip(var_with_bounds, bounds_variables):
        f.variables[var].setncattr("bounds", bndvar)

    f.close()

    return None


def parse_args():
    """parse command line arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="Input file")
    parser.add_argument(
        "-o", "--outfile", type=str, required=True, help="Output file name"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        required=False,
        help="Print detailed output",
    )
    parser.add_argument("-t", "--tagfile", action="store_true", required=False, help="")
    parser.add_argument(
        "-f",
        "--format",
        type=str,
        required=False,
        default="NETCDF3_64BIT",
        help="netcdf format for output file",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run()
