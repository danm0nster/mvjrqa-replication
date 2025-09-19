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

SIM_DATA = $(SIM_LINEAR_DATA) $(SIM_LOGISTIC_DATA) $(SIM_LORENZ_DATA) $(SIM_LORENZ96_DATA)

default: $(SIM_DATA)

# Grouped targets are useful for running make in parallel with -j flag.
# A somewhat bothersome workaround, since grouped targets are not avaiable
# on GNU make 3.81 which is default on Mac OS. On Linux or with newer
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
