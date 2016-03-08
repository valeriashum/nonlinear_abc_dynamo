void init_timeseries();
void init_timevar();
void output_timevar(const struct Field fldi,
					const double t);
void init_rate();
void output_rate(const struct Field fldi,
					const double t);
void init_isotropyT();
void output_isotropyT(const struct Field fldi, const double t);
/* Start of timeseries code */
void output_timeseries(double *w1, double *w2, double *w3, const double t);
/* End of timeseries code */
