#!/usr/bin/env python
###################################################################
#
##        Refine Atmospheric Diagnostics for CMIP6
#
####################################################################

import argparse
from glob import glob
import os
import subprocess

component_to_refine_no_redux = ["atmos", "aerosol"]
component_to_refine_redux_daily_to_monthly = ["atmos_daily"]


def find_matching_files(comp, workdir, only_cmip=True):
    """return list of files matching the pattern without full path"""
    matches = []
    if only_cmip:
        full_paths = glob(os.path.join(workdir, f"*{comp}*cmip*.nc"))
    else:
        full_paths = glob(os.path.join(workdir, f"*{comp}*.nc"))

    for f in full_paths:
        matches.append(f.split(os.path.sep)[-1])
    return matches


def fileout_name(filein):
    """build the output filename for refined"""

    parts = filein.split(".")
    # typically 2nd string has the name
    parts[1] = parts[1].replace("_cmip", "") + "_refined"
    sep = "."
    fileout = sep.join(parts)
    return fileout


def fileout_name_redux_daily_to_month(filein):
    """build the output filename for refined"""

    parts = filein.split(".")
    # typically 2nd string has the name
    parts[1] = parts[1].replace("daily", "month").replace("_cmip", "") + "_refined"
    sep = "."
    fileout = sep.join(parts)
    return fileout


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--source_dir",
        type=str,
        required=True,
        help="source directory for refine scripts",
    )
    parser.add_argument(
        "-r",
        "--refineDiagDir",
        type=str,
        required=True,
        help="output directory for refined files",
    )
    parser.add_argument(
        "-w",
        "--workdir",
        type=str,
        required=False,
        default=".",
        help="working directory",
    )
    parser.add_argument(
        "--only_cmip",
        action="store_true",
        required=False,
        help="only treat files labelled cmip",
    )
    args = parser.parse_args()

    # --- no reduction in time
    files_to_process_no_redux = []
    for comp in component_to_refine_no_redux:
        files_to_process_no_redux += find_matching_files(
            comp, args.workdir, args.only_cmip
        )

    files_to_process_no_redux.sort()

    for f in files_to_process_no_redux:
        fout = fileout_name(f)
        out = subprocess.check_call(
            [
                f"{args.source_dir}/refine_Atmos_no_redux.py",
                f,
                "-o",
                f"{args.refineDiagDir}/{fout}",
                "-v",
            ],
            shell=False,
        )

    # --- reduce from daily to monthly
    files_to_process_redux_daily_to_monthly = []
    for comp in component_to_refine_redux_daily_to_monthly:
        files_to_process_redux_daily_to_monthly += find_matching_files(
            comp, args.workdir, args.only_cmip
        )

    files_to_process_redux_daily_to_monthly.sort()

    for f in files_to_process_redux_daily_to_monthly:
        fout = fileout_name_redux_daily_to_month(f)
        out = subprocess.check_call(
            [
                f"{args.source_dir}/refine_Atmos_redux_daily_to_monthly.py",
                f,
                "-o",
                f"{args.refineDiagDir}/{fout}",
                "-v",
            ],
            shell=False,
        )
