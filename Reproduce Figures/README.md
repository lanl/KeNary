#### This folder contains the files used to reproduce the figures presented in Stricklin et al. (2024). Note: Due to storage limitations induced by Los Alamos National Laboratory, results from experiments must be downloaded from [this google drive](https://drive.google.com/drive/folders/1giX6cbTxLFEvis-Q_Y8DHdaFZ8IGcDcI?usp=sharing).

The following .R files can be used to *reproduce the figures* presented in Stricklin et al. (2024):
1. 0_SubmitFigs.R: contains the R code to reproduce the figures; 
2. KeNaryFuns.R: contains the functions used to run KeNary, reproduce figures, and simulate spectra.

The following .csv files contain *data* used in Section 4: Method Demonstration: A Forensic Paint Dataset of Stricklin et al. (2024):
1. paintcan_19.csv: chemical data corresponding to paint can number 19;
2. paintcan_34.csv: chemical data corresopnding to paint can number 34; 

The following .rData file contains the data used to *reproduce Figure 5* in Stricklin et al. (2024):
1. RMP_dat_N5_M3_fixed_trace.rData: contains the data used to make the RMP plot presented in Figure 5.

The following .csv file contains the sources of the paint data studied in each of the experiments considered in *Section 4: Method Demonstration: A Forensic Paint Dataset*:
1. source_grid.csv: source combinations considered in each experiment in Section 4.

The following .h5 files contain the data used in *Section 5: A Blind Study: A Forensic Soot Dataset* of Stricklin et al. (2024):
1. CompB_Air_1d_V_BC.h5: chemical data corresponding to soot formed under air atmosphere during detonation of composition B-3; 
2. CompB_Ar_1d_V_BC.h5: chemical data corresponding to soot formed under argon atmosphere during detonation of composition B-3; 
3. CompB_Ful_1d_V_BC.h5: chemical data corresponding to laboratory standard Fullerene soot;
4. D160707.h5: chemical data associated with a blind sample.

>[!NOTE]
> Due to storage limitations induced by Los Alamos National Laboratory, results from experiments presented in *Section 4: Method Demonstration: A Forensic Paint Dataset* must be downloaded from [this google drive](https://drive.google.com/drive/folders/1giX6cbTxLFEvis-Q_Y8DHdaFZ8IGcDcI?usp=sharing). Click on the link and request access, and M.A. Stricklin will grant you access to the data. In particular, the user must access this link to access the following files:
>
> 1. Sim_Spectra.rda: simulated spectra for all paint cans used for various experiments;
> 2. FTIRSimsResults_*.rda: results corresponding to 2-8 source scenario experiments.
