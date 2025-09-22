# Groups of files to use in targets
SIM_LINEAR_DATA = data/linear_stochastic_coupling_sweep.csv \
	data/linear_stochastic_extreme_coupling.csv \
	data/linear_stochastic_rr_sweep.csv

SIM_LOGISTIC_DATA = data/driven_coupled_logistic_maps_coupling_sweep.csv \
	data/driven_coupled_logistic_maps_extreme_coupling.csv \
	data/driven_coupled_logistic_maps_rr_sweep.csv

SIM_LORENZ_DATA = data/lorenz_harmonic_coupling_sweep.csv \
	data/lorenz_harmonic_extreme_coupling.csv \
	data/lorenz_harmonic_rr_sweep.csv

SIM_LORENZ96_DATA = data/lorenz96_harmonic_coupling_sweep.csv \
	data/lorenz96_harmonic_extreme_coupling.csv \
	data/lorenz96_harmonic_rr_sweep.csv

SIM_HIDIM_DATA = data/lorenz96_high_dim_coupling_sweep.csv \
	data/lorenz96_high_dim_extreme_coupling.csv \
	data/lorenz96_high_dim_rr_sweep.csv

SIM_MDRQA_DATA = data/linear_stochastic_system_mdrqa_comparison.csv \
	data/logistic_system_mdrqa_comparison.csv \
	data/lorenz_harmonic_mdrqa_comparison.csv \
	data/lorenz_96_harmonic_mdrqa_comparison.csv

# All of the simulated data in one variable
MODEL_DATA = $(SIM_LINEAR_DATA) $(SIM_LOGISTIC_DATA) $(SIM_LORENZ_DATA) $(SIM_LORENZ96_DATA)
SIM_DATA = $(MODEL_DATA) $(SIM_HIDIM_DATA) $(SIM_MDRQA_DATA)

PLT_MODEL_RESULTS = plots/linear_stochastic_jrr_rr.pdf \
	plots/logistic_driven_jrr_rr.pdf \
	plots/lorenz_harmonic_jrr_by_rr.pdf \
	plots/lorenz96_harmonic_jrr_rr.pdf \
	plots/linear_jrr_by_rr2_log_color.pdf \
	plots/logistic_maps_external_driver_color.pdf \
	plots/lor-osc_jrr_by_rr2_log_color.pdf \
	plots/lor96_osc_jrr_by_rr2_log_color.pdf

PLT_MDRQA = plots/comparison_mdrqa_mvjrqa_coupling.pdf

PLT_HIDIM = plots/lorenz96_high_dim_time_series_plot.pdf \
	plots/lorenz96_high_dim_harmonic_jrr_rr.pdf \
	plots/lor_hi_dim_jrr_by_rr2_log_color.pdf \
	plots/lorenz96_high_dim_extreme_coupling.pdf

ANIMATION = plots/mvjrp_animation_snapshot.pdf \
	plots/mvjrp_animation.mp4

# All plots based on simulation data
SIM_PLOTS = $(PLT_MODEL_RESULTS) $(PLT_MDRQA) $(PLT_HIDIM) $(ANIMATION)

default: data

data: $(SIM_DATA)

plots: $(SIM_PLOTS)

clean:
	rm .linear .logistic .lorenz .lorenz96 .hidim .mdrqa

# Rules for making (groups of) files from R-scripts.

# Grouped targets are useful for running make in parallel with -j flag.
# Since grouped targets are not avaiable on GNU make 3.81 which is the
# default on Mac OS, we have to use a workaround:
# Each group of files is set to depend on a hidden file beginning
# with a ".". This file then has a rule to run the required script
# and afterwards the hidden file is created or has its time stamp
# updated with the touch command.
#
# A somewhat bothersome workaround, On Linux or with newer
# make installed manually on Mac OS, we could simply use:
# $(SIM_LINEAR_DATA) &: generate_linear_stochastic_data.R
#	Rscript $<

$(SIM_LINEAR_DATA): .linear ;

.linear: generate_linear_stochastic_data.R
	Rscript $<
	@touch $@

$(SIM_LOGISTIC_DATA): .logistic ;

.logistic: generate_driven_coupled_logistic_maps_data.R
	Rscript $<
	@touch $@

$(SIM_LORENZ_DATA): .lorenz ;

.lorenz: generate_lorenz_harmonic_oscillator_data.R
	Rscript $<
	@touch $@

$(SIM_LORENZ96_DATA): .lorenz96 ;

.lorenz96: generate_lorenz96_high_dimensional_oscillator_data.R
	Rscript $<
	@touch $@

$(SIM_HIDIM_DATA): .hidim ;

.hidim: generate_lorenz96_high_dimensional_oscillator_data.R
	Rscript $<
	@touch $@

$(SIM_MDRQA_DATA): .mdrqa ;

.mdrqa: generate_data_for_comparison_to_mdrqa.R
	Rscript $<
	@touch $@

$(PLT_MODEL_RESULTS): $(MODEL_DATA) .plt_model;

.plt_model: plot_model_results.R
	Rscript $<
	@touch $@

$(PLT_MDRQA): $(SIM_MDRQA_DATA) .plt_mdrqa ;

.plt_mdrqa: plot_comparison_to_mdrqa.R
	Rscript $<
	@touch $@

$(PLT_HIDIM): $(SIM_HIDIM_DATA) .plt_hidim ;

.plt_hidim: plot_high_dimensional_lorenz96.R
	Rscript $<
	@touch $@

$(ANIMATION): .animation ;

.animation: create_mvjrp_animation.R
	Rscript $<
	@touch $@