[![CC BY 4.0][cc-by-shield]][cc-by]

# Replication data and code for _MvJRQA_ 

This repository contains R code to replicate the results in the article

> Sebastian Wallot & Dan Mønster
>
> Multivariate Joint Recurrence Quantification Analysis:
> detecting coupling between time series of different dimensionalities


## Data availability
All of the synthetic model datasets are included in this repository and scripts
to reproduce the datasets are provided (see below). The output from the
empirical data analysis is also included, and can also be reproduced if the
original data are downloaded.

The dataset for the empirical example used in the article is openly available
(Shafiei & Shadpour, 2023), but not included here. It is available through
[PhysioNet](https://physionet.org/content/eeg-eye-gaze-for-fls-tasks/1.0.0/)
(Goldberger et al., 2000), and described in an original publication
(Shafiei et al., 2023).

Note: There are errors in the file `PerformanceScores.csv` from the PhysioNet
dataset that were not corrected at the time of publication, so we provide a
corrected version here, named `PerformanceScores_fixed.csv` which is the
version used in the scripts for analyzing the empirical data.

## List of scripts and their use

The main function mvjrqa() to perform MvJRQA analysis is contained in the
file `R/mvjrqa.R`. It is used in the scripts to generate the data used in the
article. Additional functions used in the scripts are all located in files
in the directory `R`. We also rely on many packages — first and foremost
on the crqa package (Coco et al., 2021). See the section _Software and
packages used_.

### Scripts

Below is a description of the scripts included in this repository. They
are divided into three categories:

* Data generating scripts. These scripts compute the data for the
  synthetic models. You must run these scripts (or use the provided
  pre-computed data) before you can use the scripts for creating figures.
* Scripts for producing figures. These scripts use the generated data
  to produce the plots used in the figures in the article. Additionally,
  one script produces an animation available as supplementary information.
* Scripts to analyze empirical data. These scripts rely on data not 
  included here. In order to run these scripts you have to first
  download the empirical dataset (see: Shafiei & Shadpour, 2023).
  See more details below in the section _Scripts to analyze empirical data_.
  
The data generating scripts use random initial conditions, so in
order to be reproducible, the random seed is set in each script
and they run in a single thread.

#### Running the scripts individually
In order to avoid unintended side effects it is safest to run the
scripts in a new R session. One simple way of doing this is to use
`Rscript` (Unix, Linux, Mac OS) or `Rscript.exe` (Windows). For
example to run the first data generation script below, open a terminal
change to the directory where the script is located and then run:

```
Rscript generate_linear_stochastic_data.R
```

#### Running the scripts with make
A Makefile has been included to automate running the required scripts
for generating data from simulated models and for creating the plots 
used in the figures as well as the movie. This should work on Mac OS,
Linux and other Unix systems.

All the simulation data and plots are already included in the
sub-directories `data` and `plots` respectively. In order to replicate these
files from the scripts provided you can simply type `make data` and 
`make plots`. Note that _this will overwrite the original files_. 

To speed up the execution, you can process several scripts in parallel,
by providing the j-flag to make along with the number of parallel processes
to use. For example, using 6 CPU cores you can write:

```
make -j6 data
```

and

```
make -j6 plots
```

#### Data generating scripts

The list below shows the scripts that generate data and the output files
generated. Running the scripts again will overwrite the original files.

* generate_linear_stochastic_data.R
  - data/linear_stochastic_coupling_sweep.csv
  - data/linear_stochastic_extreme_coupling.csv
  - data/linear_stochastic_rr_sweep.csv
* generate_driven_coupled_logistic_maps_data.R
  - data/driven_coupled_logistic_maps_coupling_sweep.csv
  - data/driven_coupled_logistic_maps_extreme_coupling.csv
  - data/driven_coupled_logistic_maps_rr_sweep.csv
* generate_lorenz_harmonic_oscillator_data.R
  - data/lorenz_harmonic_coupling_sweep.csv
  - data/lorenz_harmonic_extreme_coupling.csv
  - data/lorenz_harmonic_rr_sweep.csv
* generate_lorenz96_harmonic_oscillator_data.R
  - data/lorenz96_harmonic_coupling_sweep.csv
  - data/lorenz96_harmonic_extreme_coupling.csv
  - data/lorenz96_harmonic_rr_sweep.csv
* generate_lorenz96_high_dimensional_oscillator_data.R
  - data/lorenz96_high_dim_coupling_sweep.csv
  - data/lorenz96_high_dim_extreme_coupling.csv
  - data/lorenz96_high_dim_rr_sweep.csv
* generate_data_for_comparison_to_mdrqa.R
  - data/linear_stochastic_system_mdrqa_comparison.csv
  - data/logistic_system_mdrqa_comparison.csv
  - data/lorenz_96_harmonic_mdrqa_comparison.csv
  - data/lorenz_harmonic_mdrqa_comparison.csv

#### Scripts for producing figures

Run these scripts individually using Rscript or use `make plots` to 
run the scripts in series or `make -j6 plots` to run them in parallel.

* plot_model_time_series.R
* plot_model_results.R
* plot_comparison_to_mdrqa.R
* plot_extreme_coupling.R
* plot_high_dimensional_lorenz96.R
* plot_phase_space_coupling.R
* plot_noise_effect.R
* plot_conceptual_jrci.R
* create_mvjrp_animation.R

#### Scripts to analyze empirical data
There are three scripts to analyze the empirical dataset and they must be
run in the following order:

* 01_empirical_mvjrqa_analysis.R
* 02_empirical_merge_files.R
* 03_empirical_mvjrqa_regression.R
* 04_empirical_plots.R

The first script requires more resources to run, since it loops over all the
observations and performs MvJRQA for EEG and 2D eye tracking as well as EEG
and 3D eye tracking data. This script uses parallel processing and you have
to set the number of worker processes in the script if you want to change the
default to tailor it to your own computing environment (see notes in the 
script header). **Note:** This script can take a long time to run and requires
a sufficient amount of memory (see notes in the script header).

The second script merges the data created by the first script into the file
`results/mvjrqa.csv`. 

The third script reads in the merged data and performs a regression analysis.
This script also plots the model results and produces tables of the results.

The fourth scripts creates plots of the joint recurrence coupling indicator,
example time series for a single experimental trial, and some diagnostic
plots for the supplementary information.

## Software and packages used
The scripts use R (R Core Team, 2025) and the following packages:

* animation: Xie, 2013; Xie et al., 2021.
* cowplot: Wilke, 2024.
* crqa: Coco et al., 2021.
* data.table: Barret et al., 2025.
* deSolve: Soetaert et al., 2010.
* dplyr: Wickham et al., 2021.
* edf: Henelius, 2016.
* effectsize: Mattan et al., 2020.
* gridExtra: Auguie, 2017.
* ggplot2: Wickham, 2016.
* ggrepel: Slowikowski, 2024.
* latex2exp: Meschiari, 2022.
* lme4: Bates et al., 2015.
* lmerTest: Kuznetsova et al., 2017.
* Matrix: Bates et al., 2025.
* modelbased: Makowski et al, 2025.
* parallel: R Core Team, 2025.
* patchwork: Pedersen, 2024.
* performance: Lüdecke et al., 2021a.
* purrr: Wickham & Henry, 2025.
* readr: Wickham et al., 2024a.
* rgl: Murdoch & Adler, 2025.
* see: Lüdecke et al., 2021b.
* texreg: Leifeld, 2023.
* tidyr: Wickham et al., 2024b.

In order to reproduce the animation video, you also need to have
[FFmpeg](https://www.ffmpeg.org/) installed on your system.

## References

Auguie B (2017). _gridExtra: Miscellaneous Functions for "Grid" Graphics_.
  R package version 2.3, <https://CRAN.R-project.org/package=gridExtra>
  
Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T,
  Schwendinger B, Krylov I (2025).
  _data.table: Extension of data.frame_. R package version 1.17.0,
  <https://CRAN.R-project.org/package=data.table>.
  
Bates D, Maechler M, Bolker B, Walker S (2015).
  _Fitting Linear Mixed-Effects Models Using lme4._
  Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
  
Bates D, Maechler M, Jagan M (2025). _Matrix: Sparse and Dense Matrix Classes
  and Methods_. R package version 1.7-2,
  <https://CRAN.R-project.org/package=Matrix>.
  
Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect
  Size Indices and Standardized Parameters. Journal of Open Source Software,
  5(56), 2815. <https://doi.org/10.21105/joss.02815>.

Coco MI, Mønster D, Leonardi G, Dale R, Wallot S (2021). _Unidimensional and
  Multidimensional Methods for Recurrence Quantification Analysis with crqa._
  R Journal, *13*(1). doi:10.32614/RJ-2021-062
  <https://doi.org/10.32614/RJ-2021-062>.
  
Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark,
  R., ... & Stanley, H. E. (2000). _PhysioBank, PhysioToolkit, and PhysioNet:
  Components of a new research resource for complex physiologic signals._
  Circulation [Online]. 101 (23), pp. e215–e220. RRID:SCR_007345.
  
Henelius A (2016). _edf: Read Data from European Data Format (EDF and EDF+)
  Files_. R package version 1.0.0, <https://CRAN.R-project.org/package=edf>.

Kuznetsova A, Brockhoff PB, Christensen RHB (2017).
  _lmerTest Package: Tests in Linear Mixed Effects Models._
   Journal of Statistical Software, *82*(13), 1-26. doi:10.18637/jss.v082.i13
  <https://doi.org/10.18637/jss.v082.i13>.

Leifeld, P. (2013).
  _texreg: Conversion of Statistical Model Output in R to LATEX and HTML Tables._
  Journal of Statistical Software, 55, 1–24.
  <https://doi.org/10.18637/jss.v055.i08>

Lüdecke D,  Ben-Shachar M S, Patil I, Waggoner P, Makowski P, (2021a).
  _performance: An R Package for Assessment, Comparison and Testing of
  Statistical Models._
  Journal of Open Source Software, *6*(60), 3139.
  <https://doi.org/10.21105/joss.03139>.

Lüdecke, D., Patil, I., Ben-Shachar, M. S., Wiernik, B. M., Waggoner, P.,
  & Makowski, D. (2021b). 
  _see: An R Package for Visualizing Statistical Models._
  Journal of Open Source Software, 6(64), 3393.
  <https://doi.org/10.21105/joss.03393>

Makowski, D., Ben-Shachar, M. S., Wiernik, B. M., Patil, I., Thériault,
  R., & Lüdecke, D. (2025). _modelbased: An R package to make the most out of
  your statistical models through marginal means, marginal effects, and
  model predictions._
  Journal of Open Source Software, 10(109), 7969.
  <https://doi.org/10.21105/joss.07969>

Meschiari S (2022). _latex2exp: Use LaTeX Expressions in Plots_.
  R package version 0.9.6, <https://CRAN.R-project.org/package=latex2exp>.
  
Murdoch D, Adler D (2025). _rgl: 3D Visualization Using OpenGL_.
  R package version 1.3.18, <https://CRAN.R-project.org/package=rgl>
  
Pedersen T (2024). _patchwork: The Composer of Plots_. R package version 1.3.0,
  <https://CRAN.R-project.org/package=patchwork>.
  
R Core Team (2025). _R: A Language and Environment for Statistical Computing_.
  R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.
  
Shafiei, S. B., & Shadpour, S. (2023).
  _Integration of Electroencephalogram and Eye-Gaze Datasets for Performance
  Evaluation in Fundamentals of Laparoscopic Surgery (FLS) Tasks_
  (version 1.0.0). PhysioNet. RRID:SCR_007345.
  https://doi.org/10.13026/kyjw-p786.
  
Somayeh B. Shafiei, Saeed Shadpour, Xavier Intes, Rahul Rahul,
  Mehdi Seilanian Toussi, Ambreen Shafqat, _Performance and Learning Rate
  Prediction Models Development in FLS and RAS Surgical Tasks Using
  Electroencephalogram and Eye Gaze Data and Machine Learning,_
  Surgical Endoscopy, 2023.

Slowikowski K (2024). 
  _ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'_.
  R package version 0.9.6, <https://CRAN.R-project.org/package=ggrepel>.
  
Soetaert Karline, Thomas Petzoldt, R. Woodrow Setzer (2010).
  Solving Differential Equations in R: Package deSolve.
  Journal of Statistical Software, 33(9), 1--25. doi:10.18637/jss.v033.i09

Wickham, Hadley. _ggplot2: Elegant Graphics for Data Analysis_
  Springer-Verlag New York, 2016.

Wickham, Hadley, Romain François, Lionel Henry and Kirill Müller (2021). _dplyr:
  A Grammar of Data Manipulation._  https://CRAN.R-project.org/package=dplyr
  
Wickham H, Bryan J (2025). _readxl: Read Excel Files_. R package version 1.4.5,
  <https://CRAN.R-project.org/package=readxl>.

Wickham H, Hester J, Bryan J (2024a). _readr: Read Rectangular Text Data_.
  R package version 2.1.5, <https://CRAN.R-project.org/package=readr>.
  
Wickham H, Vaughan D, Girlich M (2024b). _tidyr: Tidy Messy Data_.
  R package version 1.3.1, <https://CRAN.R-project.org/package=tidyr>.
  
Wickham H, Henry L (2025). _purrr: Functional Programming Tools_.
  R package version 1.0.4, <https://CRAN.R-project.org/package=purrr>.

Wilke C (2024). _cowplot: Streamlined Plot Theme and Plot Annotations for
  'ggplot2'_. R package version 1.1.3,
  <https://CRAN.R-project.org/package=cowplot>.
  
Xie Yihui (2013). _animation: An R Package for Creating Animations and 
  Demonstrating Statistical Methods._ Journal of Statistical Software,
  53(1), 1-27. URL https://doi.org/10.18637/jss.v053.i01.

Xie Y, Mueller C, Yu L, Zhu W (2021). _animation: A
  Gallery of Animations in Statistics and Utilities to Create Animations._
  R package version 2.7, <https://yihui.org/animation>.

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
