# Script for MultiTaper Method-Singular Value Decomposition (MTM-SVD) with Monte Carlo test in python


## outlook of the next version

The shortage of the version is being too slow.
Making the function running with numba may help

---

This script is a modified version of the Python function developed by
Mathilde Jutras at McGill University, Canada[1]. 
You can find the original Python code here: 
https://github.com/mathildejutras/mtm-svd-python

This script was adapted by Yitao Liu at Nanjing University of Information Science & Technology, China
Copyright (C) 2021, Yitao Liu
and is available under the GNU General Public License v3.0

The script may be used, copied, or redistributed as long as it is cited as follow:
[]

This software may be used, copied, or redistributed as long as it is not 
sold and that this copyright notice is reproduced on each copy made. 
This routine is provided as is without any express or implied warranties.

Questions or comments to:
Yitao Liu, liuyitao97@outlook.com

---

The script is structured as follows:

In the main script is found in mtm-svd-python.py
In the first section, the user can load the data,
assuming the outputs are stored in a netcdf format.
In the secton section, functions are called to calculate the spectrum
The user will then be asked for which frequencies he wants to plot 
the spatial patterns associated with the variability.
In the third section, the spatial patterns are plotted and saved

The main script is contained in mtm-svd-python.py, and the required functions can be found in mtm_functions.py.

---

Python Package needed:
- numpy
- scipy
- xarray (read the netcdf file)
- matplotlib (not necessary, just for plotting)

You can install the needed Python packages by conda,with the command below
```
conda install -c conda-forge numpy scipy xarray matplotlib
```

[1] Mathilde Jutras. (2020, July 6). mathildejutras/mtm-svd-python: v1.0.0-alpha (Version v1.0.0). Zenodo. http://doi.org/10.5281/zenodo.3932319

