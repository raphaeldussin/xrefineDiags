#!/usr/bin/env python

# This script contains the refineDiags that reduce data in time
# from daily to monthly values (e.g. min/ax over the period)
# It can accept any file and will only compute refineDiags in fields
# are present.

import argparse
import cftime
import os
import netCDF4 as nc
import numpy as np
import xarray as xr

xr.set_options(keep_attrs=True)

CMOR_MISSING_VALUE = 1.0e20
extra_time_variables = ["time_bnds", "average_T1", "average_T2", "average_DT"]
do_not_encode_vars = ["nv", "grid_xt", "grid_yt", "time", "lev"]
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

    ds = xr.open_dataset(args.infile, use_cftime=True)
    # make a copy of undecoded time
    time_cop = xr.open_dataset(args.infile, decode_cf=False)["time"]

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

    # we haven't created any new variables yet
    new_vars_output = False

    # -- create time variables if they don't exist in refined dataset
    if "time" not in list(refined.variables):
        refined = populate_time_vars(refined, ds, time_cop)
        # creating time variables do not qualify as new output
        new_vars_output = False

    # -- compute min/max from daily to monthly
    if "t_ref_max" in list(ds.variables):
        new_vars_output = True
        refined["t_ref_max"] = max_over_month(ds["t_ref_max"], ds, refined)
        refined["t_ref_max"].encoding["_FillValue"] = CMOR_MISSING_VALUE
    if "t_ref_min" in list(ds.variables):
        new_vars_output = True
        refined["t_ref_min"] = min_over_month(ds["t_ref_min"], ds, refined)
        refined["t_ref_min"].encoding["_FillValue"] = CMOR_MISSING_VALUE

    if "u_ref_max" in list(ds.variables):
        new_vars_output = True
        refined["u_ref_max"] = max_over_month(ds["u_ref_max"], ds, refined)
        refined["u_ref_max"].encoding["_FillValue"] = CMOR_MISSING_VALUE
    if "u_ref_min" in list(ds.variables):
        new_vars_output = True
        refined["u_ref_min"] = min_over_month(ds["u_ref_min"], ds, refined)
        refined["u_ref_min"].encoding["_FillValue"] = CMOR_MISSING_VALUE

    if "rh_ref_max" in list(ds.variables):
        new_vars_output = True
        refined["rh_ref_max"] = max_over_month(ds["rh_ref_max"], ds, refined)
        refined["rh_ref_max"].encoding["_FillValue"] = CMOR_MISSING_VALUE
    if "rh_ref_min" in list(ds.variables):
        new_vars_output = True
        refined["rh_ref_min"] = min_over_month(ds["rh_ref_min"], ds, refined)
        refined["rh_ref_min"].encoding["_FillValue"] = CMOR_MISSING_VALUE

    # --- write dataset to file
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


def cftime_to_float(time_cf, time_cop):
    """ convert the cftime representation into float using
    time_cop as a template for units, epoch,... """

    # get units and calendar
    units = time_cop.attrs["units"]

    if "calendar_type" in list(time_cop.attrs):
        cal = time_cop.attrs["calendar_type"].lower()
    elif "calendar" in list(time_cop.attrs):
        cal = time_cop.attrs["calendar"].lower()
    else:
        raise ValueError("No calendar attributes for time variable")

    time_noncf = cftime.date2num(time_cf, units, cal)

    return time_noncf

def populate_time_vars(refined, ds, time_cop):

    realtime = xr.decode_cf(ds)["time"]

    # compute monthly values
    time_monthly = (
        ds["time"].groupby(realtime.dt.month).mean(dim="time").rename({"month": "time"})
    )

    time_noncf = cftime_to_float(time_monthly, time_cop)

    refined["time"] = xr.DataArray(
        data=time_noncf, dims=("time"), attrs=time_cop.attrs
    )

    average_DT = (
        ds["average_DT"]
        .groupby(realtime.dt.month)
        .sum(dim="time")
        .rename({"month": "time"})
    )
    refined["average_DT"] = xr.DataArray(data=average_DT.values.astype("f8"), dims=("time"))
    refined["average_DT"].attrs = ds["average_DT"].attrs.copy()
    refined["average_DT"].attrs["units"] = "days"
    refined["average_DT"].encoding["_FillValue"] = CMOR_MISSING_VALUE

    average_T1 = (
        ds["average_T1"]
        .groupby(realtime.dt.month)
        .min(dim="time")
        .rename({"month": "time"})
    )

    average_T1_noncf = cftime_to_float(average_T1, time_cop)
    refined["average_T1"] = xr.DataArray(data=average_T1_noncf.astype("f8"), dims=("time"))
    refined["average_T1"].attrs = ds["average_T1"].attrs.copy()
    refined["average_T1"].attrs["units"] = time_cop.attrs["units"]
    refined["average_T1"].encoding["_FillValue"] = CMOR_MISSING_VALUE

    average_T2 = (
        ds["average_T2"]
        .groupby(realtime.dt.month)
        .max(dim="time")
        .rename({"month": "time"})
    )

    average_T2_noncf = cftime_to_float(average_T2, time_cop)
    refined["average_T2"] = xr.DataArray(data=average_T2_noncf.astype("f8"), dims=("time"))
    refined["average_T2"].attrs = ds["average_T2"].attrs.copy()
    refined["average_T2"].attrs["units"] = time_cop.attrs["units"]
    refined["average_T2"].encoding["_FillValue"] = CMOR_MISSING_VALUE

    time_bnds_noncf = np.stack([average_T1_noncf.astype("f8"), average_T2_noncf.astype("f8")], axis=-1)
    refined["time_bnds"] = xr.DataArray(data=time_bnds_noncf, dims=("time", "nv"))
    refined["time_bnds"].attrs = ds["time_bnds"].attrs.copy()
    refined["time_bnds"].attrs["units"] = time_cop.attrs["units"]
    refined["time_bnds"].encoding["_FillValue"] = CMOR_MISSING_VALUE

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
