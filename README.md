# YT-Simulation: Data Visualization

This repository provides tools for visualizing astronomical datasets using the `yt` package, specifically optimized for JWST filter integration.

## Quick Start:
To simulate these datasets without errors, you will need three things: **filter profiles**, a **dataset**, and the **coding notebooks**.

### 1. JWST Filter Profile
This project requires a JWST filter file.
* **Custom Filters:** Download official profiles from the [NIRCam Filters website](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters).
* **Included Example:** If you don't know what filter you want to try with this code, I have provided the [F200W filter txt file](./F200W_filter.txt) in this repo as an example.

### 2. Research Dataset
This code assumes you have an object dataset. 
* If you do not have your own simulation data and want to make sure that your dataset format matches the format used in this code, you can download my sample dataset here: [Raw data v1.0](https://github.com/Patsrnpt/YT-Simulation/releases/tag/v1.0).

Because the raw data folder size is very large, it exceeds GitHub's standard upload limits. To resolve this, I have compressed the data into multiple 1GB .zip files. In order to use these files, you should:

1. Download all split parts (**output_part_aa.zip** through **output_part_ae.zip**) from the release link above.
2. Combine them into a single archive:
   * **On Mac/Linux:** Use the terminal command `cat output_part_*.zip > output_00273.zip`
   * **On Windows:** Use the command `copy /b output_part_*.zip output_00273.zip`
3. Keeping these datasets in a single folder will ensure the analysis scripts can locate the data correctly.

### 3. Coding Files

This repository includes 2 different parts to process the data:

* **[Library Folder](./library)**: This folder includes all .py files used for visualization. To access this code, navigate to the library directory. It contains 6 specialized files that handle everything from math to visualization:
    * **`calculate_quantities.py`**: Handles core mathematical transformations (flux, magnitude, and surface brightness).
    * **`continuum_grid.py`**: Uses **PyNeb** to generate physics-based continuum grids (bound-free, free-free, two-photon).
    * **`filter_tools.py`**: Contains utilities for generating and managing wavelength filter bins.
    * **`prepare_data.py`**: Manages dataset loading, star particle centering, and coordinate transformations.
    * **`visualization.py`**: Plotting tools for creating Projection plots, Phase plots, and Spectra.
    * **`yt_fields.py`**: Integrates the calculated quantities into the `yt` package as derived fields.

* **'visualization_example.ipynb'**: This Jupyter notebook provides a command of how to simulate a projection plot, a phase plot, and a spectrum of your dataset using the functions provided in the library folder.

---

### How to Import
To use these functions in your notebook, you must import them from the library folder:

```python
# import all related functions from the library folder
from library.calculate_quantities import *
from library.continuum_grid import *
from library.filter_tools import *
from library.prepare_data import *
from library.visualization import *
from library.yt_fields import *
```

If you have any questions, feel free to contact me via this email: [sphoom22@terpmail.umd.edu](mailto:sphoom22@terpmail.umd.edu)
