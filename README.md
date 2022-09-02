# xrefineDiags: a xarray-based implementation of the refineDiags

## Install:

Clone the repo!

## Usage:

the individual python scripts are executable and can take a single netcdf file as argument and produce refined output file.
For convenience, several meta-scripts and wrappers are available to batch process all files.
For example, `refineDiag_{atmos,ocean}.csh` and `refineDiag_{atmos,ocean}.py` allow to process all atmospheric files.

```bash
refineDiag_atmos.csh PATH_TO_REFINEDIAGS OUTPUT_DIRECTORY
```

```python
refineDiag_atmos.py -s PATH_TO_REFINEDIAGS -r OUTPUT_DIRECTORY --only_cmip
```

To call this from inside the xml, add the following call to the FRE wrapper :

```xml
<refineDiag script="PATH_TO_INSTALL/xrefineDiags/refineDiag_atmos_FRE.csh"/>
```
