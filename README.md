# LinaForma
<p align="center">
  <!-- DOI -->
<a href="https://zenodo.org/doi/10.5281/zenodo.11110440">
  <img src="https://zenodo.org/badge/790819288.svg" alt="DOI">
</a>
  <!-- License -->
  <a href="https://www.gnu.org/licenses/gpl-3.0">
    <img src="https://img.shields.io/badge/License-GPLv3-blue.svg" />
  </a>
</p>

 <p align="center">
<img src="https://github.com/TMackay-Champion/LinaForma/blob/05e58a21e651066dc0452beaa799e8eab52530d0/images/logo_heatmap.jpg", width="60%">
</p>


LinaForma is a series of MATLAB® scripts for calculating the optimal pressure-temperature (P-T) conditions experienced by a rock using a grid-search inversion, accompanied by bootstrap re-sampling to quantify the solution uncertainty and sensitivity to the input variables. LinaForma calculates the difference ("misfit") between measurements (e.g., Xalm) and forward models computed for a range of P-T points in third-party software such as [THERIAK-DOMINO](https://titan.minpet.unibas.ch/minpet/theriak/prog11032020/), [Perple_X](https://www.perplex.ethz.ch/), and [MAGEMin](https://github.com/ComputationalThermodynamics/MAGEMin). The P-T point with the lowest misfit value defines the “best-fit” solution. A suite of tools is also provided for plotting forward model data, performing tasks such as Principal Component Analysis, and automating relevant processes in THERIAK-DOMINO. 

LinaForma requires no prior computer programming knowledge and a step-by-step walkthrough is provided.


Key Features
-------
There are three key packages:
1) **LinaForma_inversion**: this package performs a grid-search inversion to quantify which P-T conditions best fit the rock of interest. Bootstrap re-sampling provides an assessment of the uncertainty associated with this P-T estimate, and the sensitivity of the estimate to uncertainty in the input variables.
2) **TheriakDomino_tools**: this package provides useful tools for creating bulk-composition THERIN files (E1 script), collating DOMINO output into a single csv. file (E2), automating THERIAK (E3), plotting data in heatmaps, contour plots etc. (E4), and calculating the correlation coefficients between the different variables, temperature and pressure (E5).
3) **PCA_phase_composition**: this package allows the user to perform Principal Component Analysis and K-means clustering on mineral composition data.


Citation
--------
If you find this package useful, please do consider citing it using the following:

```console
Mackay-Champion, T., & Cawood, I. (2024). LinaForma (Version 1.0.) [Computer software]. https://doi.org/10.5281/zenodo.11110441
```


Contributing to LinaForma
----------------------------
Contributions to LinaForma are welcomed. Whether you have identified a bug or would like to request a new feature or enhancement, please reach out either directly or via the GitHub Issues panel to discuss the proposed changes.

Links to projects that have made use of LinaForma are also most welcome.


Contact
-------
Any comments/questions can be directed to:
* **Tobermory Mackay-Champion** - tmackaychampion@gmail.com
* **Ian Cawood** - ipcawood@gmail.com

License
-------
This package is written and maintained by Tobermory Mackay-Champion and Ian Cawood. It is distributed under the GPLv3 License. Please see the [LICENSE](LICENSE) file for a complete description of the rights and freedoms that this provides the user.
