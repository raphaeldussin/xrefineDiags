#!/usr/bin/env python

# This script contains the refineDiags that produce data at the same
# frequency as the input data (no reduction) such as surface albedo,
# masking fields,...
# It can accept any file and will only compute refineDiags in fields
# are present.

import argparse
import os
import xarray as xr


CMOR_MISSING_VALUE = 1.0e20
extra_time_variables = ["time_bnds", "average_T1", "average_T2", "average_DT"]
do_not_encode_vars = ["nv", "grid_xt", "grid_yt", "time"]
albedo_shortname = "albs"
albedo_metadata = dict(
    long_name="Surface Albedo", units="1.0", standard_name="surface_albedo"
)
pkgname = "xrefineDiags"
scriptname = "refine_Atmos_no_redux.py"


def run():
    """run the refineDiags"""

    # --- parse command line arguemnts
    args = parse_args()
    verbose = args.verbose

    # --- open file
    if verbose:
        print(f"{pkgname}/{scriptname}: Opening input file {args.infile}")

    ds = xr.open_dataset(args.infile, decode_cf=False)

    # --- create output dataset
    if os.path.exists(args.outfile):
        if verbose:
            print(f"{pkgname}/{scriptname}: Opening existing file {args.outfile}")

        out = xr.open_dataset(args.outfile, decode_cf=False)
    else:
        if verbose:
            print(f"{pkgname}/{scriptname}: Creating new dataset")

        out = xr.Dataset()

    # --- surface albedo
    albedo_input_vars = set([args.shortwave_down, args.shortwave_up])

    if albedo_input_vars.issubset(set(ds.variables)):
        if verbose:
            print(f"{pkgname}/{scriptname}: compute surface albedo")

        out[albedo_shortname] = compute_albedo(
            ds, swdown=args.shortwave_down, swup=args.shortwave_up
        )
    else:
        if verbose:
            print(
                f"{pkgname}/{scriptname}: surface albedo NOT computed, missing input variables"
            )

    # --- add extra time variables
    for var in extra_time_variables:
        if var in list(ds.variables):
            out[var] = ds[var]

    # --- write dataset to file
    n_vars_output = len(list(out.variables))
    if n_vars_output > 0:
        out.load()
        encoding = set_netcdf_encoding(out)
        if verbose:
            print(
                f"{pkgname}/{scriptname}: writting variables {list(out.variables)} into refined file {args.outfile} "
            )

        out.to_netcdf(args.outfile, format=args.format, encoding=encoding)
    else:
        if verbose:
            print(
                f"{pkgname}/{scriptname}: no variables created, not writting refined file"
            )


def compute_albedo(ds, swdown="rsds", swup="rsus"):
    """compute surface albedo from upwelling and downwelling shortwave radiation rsus/rsds"""

    albedo = xr.where(ds[swdown] > 0.0, ds[swup] / ds[swdown], CMOR_MISSING_VALUE)

    # copy attributes from input field
    albedo.attrs = ds[swdown].attrs.copy()
    # remove interp_method for fregrid
    if "interp_method" in albedo.attrs:
        albedo.attrs.pop("interp_method")
    # update name and units
    albedo.attrs.update(albedo_metadata)
    # add comment on how albedo was computed
    albedo.attrs.update({"comment": f"{swup}/{swdown}"})

    return albedo


def set_netcdf_encoding(ds):
    """set preferred options for netcdf encoding"""

    all_vars = list(ds.variables)
    encoding = {}

    for var in do_not_encode_vars:
        if var in all_vars:
            encoding.update({var: dict(_FillValue=None)})

    return encoding


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
        "-d",
        "--shortwave_down",
        type=str,
        required=False,
        default="rsds",
        help="name of downwards shortwave radiation",
    )
    parser.add_argument(
        "-u",
        "--shortwave_up",
        type=str,
        required=False,
        default="rsus",
        help="name of upwards shortwave radiation",
    )
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
