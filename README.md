# xrefineDiags: a xarray-based implementation of the refineDiags

## Usage:

the individual python scripts are executable and can take a single netcdf file as argument and produce refined output file.
For convenience, several meta-scripts and wrappers are available to batch process all files.
For example, `refineDiag_{atmos,ocean}.csh` and `refineDiag_{atmos,ocean}.py` allow to process all atmosphere or ocean files.

```bash
refineDiag_atmos.csh PATH_TO_REFINEDIAGS OUTPUT_DIRECTORY
```

```python
refineDiag_atmos.py -s PATH_TO_REFINEDIAGS -r OUTPUT_DIRECTORY --only_cmip
```

## In the FRE workflow:

To call this from your experiments' xml, add the following call in the `<postProcess>` block of your experiment's xml:

```xml
<refineDiag script="${includeDir}/refineDiag/refineDiag_atmos_FRE.csh"/>
```

and copy the `refineDiag_atmos_FRE.csh` script into the xml repo under `awg_include/refineDiag` or `include/refineDiag` depending on your setup.

## As a standalone tool

Clone the repo! and execute the python scripts individually, use `-h` for the usage. Typically requires an output file after `-o` and an input file.
