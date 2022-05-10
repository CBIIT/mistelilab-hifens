# Fig. 6 Scripts Instructions

The scripts in this folder are meant to be run in this order:

-   `01_data_prep_Fig6.Rmd`: This script reads the single cell data obtained from the Columbus image analysis step from the `Object_Data` folder. In addition, it reads reformatting metadata (Including the position and identity of each siRNA oligo on the imaging plates. Finally, it reformats this information it in a way that is compatible with `cellHTS2` analysis of the screen, and saves it in the `hts2_input` folder.

-   `02_hts2_analysis_Fig6.Rmd`: This scripts read the data produced by the scripts above and produces the classical `cellHTS2` outputs in the `hts2_output` folder for each of the cellular features measured by Columbus.

-   `03_plotting_Fig6_B_F.Rmd`: The final scripts read the results of the `cellHTS2` analysis and saves the figures of the manuscripts in the `output` folder. The final report is in the XYZ file.
