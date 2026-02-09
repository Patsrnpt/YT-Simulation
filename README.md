# YT-Simulation: JWST Data Visualization

This repository provides tools for visualizing astronomical datasets using the `yt` package, specifically optimized for JWST filter integration.

## Quick Start:
To simulate these datasets without errors, you will need two things: a **filter profile** and a **dataset**.

### 1. JWST Filter Profiles
This project requires a JWST filter file. 
* **Custom Filters:** Download official profiles from the [NIRCam Filters website](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters).
* **Included Example:** If you are testing the code for the first time, I have provided the [F200W filter text file](./F200W_filter.txt) in this repo.

### 2. Research Dataset
This code assumes you have an object dataset. 
* If you do not have one, you can download my sample dataset here: []().
* Place the dataset in the `/data` directory to ensure the paths in the script resolve correctly.

If you have any questions, feel free to contact me via this email: [sphoom22@terpmail.umd.edu](mailto:sphoom22@terpmail.umd.edu)
