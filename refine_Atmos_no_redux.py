#!/usr/bin/env python

# This script contains the refineDiags that produce data at the same
# frequency as the input data (no reduction) such as surface albedo,
# masking fields,...
# It can accept any file and will only compute refineDiags in fields
# are present.

import argparse
import os
import netCDF4 as nc
import xarray as xr


CMOR_MISSING_VALUE = 1.0e20
extra_time_variables = ["time_bnds", "average_T1", "average_T2", "average_DT"]
do_not_encode_vars = ["nv", "grid_xt", "grid_yt", "time"]
grid_vars = ["grid_xt", "grid_yt"]
unaccepted_variables_for_masking = ["cll", "clm", "clh"]
albedo_shortname = "albs"
surface_pressure_shortname = "ps"
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

    # --- mask variables with surface pressure
    if surface_pressure_shortname in list(ds.variables):
        vars_to_process = list(ds.variables)
        vars_to_process.remove(surface_pressure_shortname)
        for var in vars_to_process:
            plev = pressure_coordinate(ds, var, verbose=verbose)
            if plev is not None:
                if var != plev.name:
                    varout = var.replace("_unmsk", "")
                    out[varout] = mask_surface_pressure(
                        ds, var, plev, ps=surface_pressure_shortname
                    )
                    out[plev.name].attrs = ds[plev.name].attrs.copy()

    # --- write dataset to file
    new_vars_output = len(list(out.variables)) > 0

    if verbose and new_vars_output:
        print(
            f"{pkgname}/{scriptname}: writting variables {list(out.variables)} into refined file {args.outfile} "
        )
    elif verbose and not new_vars_output:
        print(
            f"{pkgname}/{scriptname}: no variables created, not writting refined file"
        )

    if new_vars_output:
        write_dataset(out, ds, args)


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


def mask_surface_pressure(ds, var, pressure_dim, ps="ps"):
    """mask data with pressure larger than surface pressure"""

    plev_extended, _ = xr.broadcast(pressure_dim, ds[var])
    ps_extended, _ = xr.broadcast(ds[surface_pressure_shortname], ds[var])
    masked = xr.where(plev_extended > ps_extended, CMOR_MISSING_VALUE, ds[var])
    masked.attrs = ds[var].attrs.copy()
    masked = masked.transpose(*ds[var].dims)

    return masked


def pressure_coordinate(ds, varname, verbose=False):
    """check if dataArray has pressure coordinate fitting requirements
    and return it"""

    pressure_coord = None

    for dim in list(ds[varname].dims):
        if dim in list(ds.variables):  # dim needs to have values in file
            if ds[dim].attrs["long_name"] == "pressure":
                pressure_coord = ds[dim]
            elif ("coordinates" in ds.attrs) and (ds[dim].attrs["units"] == "Pa"):
                pressure_coord = ds[dim]

    # some variables need not to be masked
    if varname in unaccepted_variables_for_masking:
        pressure_coord = None

    if verbose:
        if pressure_coord is not None:
            print(f"{varname} has pressure coords {pressure_coord.name}")
        else:
            print(f"{varname} has no pressure coords")

    return pressure_coord


def write_dataset(ds, template, args):
    """prepare the dataset and dump into netcdf file"""

    ds.attrs = template.attrs.copy()  # copy global attributes

    # --- add proper grid attrs
    for var in grid_vars:
        if var in list(template.variables):
            ds[var] = template[var]
            ds[var].attrs = template[var].attrs.copy()

    # --- add extra time variables
    for var in extra_time_variables:
        if var in list(template.variables):
            ds[var] = template[var]
            ds[var].attrs = template[var].attrs.copy()

    encoding = set_netcdf_encoding(ds)
    ds.to_netcdf(
        args.outfile, format=args.format, encoding=encoding, unlimited_dims="time"
    )
    post_write(args.outfile, ds)

    return None


def set_netcdf_encoding(ds):
    """set preferred options for netcdf encoding"""

    all_vars = list(ds.variables)
    encoding = {}

    for var in do_not_encode_vars:
        if var in all_vars:
            encoding.update({var: dict(_FillValue=None)})

    return encoding


def post_write(filename, ds):
    """fix a posteriori attributes that xarray.to_netcdf
    did not do properly using low level netcdf lib"""

    f = nc.Dataset(filename, "a")
    f.variables["time_bnds"].setncattr("units", ds["time_bnds"].attrs["units"])
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
