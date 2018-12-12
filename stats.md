</br>
#### Interpretation

- 'Mean' and 'SD' (standard deviation) are calculated for each parameter using all post-burnin samples for either 1 or 2 chains. 

- 'ESS' stands for effective sample size and is calculated using the `effectiveSize` function of [coda](https://cran.r-project.org/web/packages/coda/index.html). The results will therefore differ from the values calculated with Phylobayes using `tracecomp`. If 2 trace files are analysed, effective sample sizes will be summed between those. *[recommended value >=100]*

- 'Geweke' stands for Geweke's convergence diagnostic and is calculated with the `geweke.diag` function of [coda](https://cran.r-project.org/web/packages/coda/index.html). Briefly, it tests for equality of the means between the first 10% and last 50% of the chain. If the chain has reached stationarity and the burnin is adequate, these means should be equal for all parameters. This can be useful in assessing if a chain has reached stationarity, and also to determine adequate burnin sizes. The value returned is a z-score. *[recommended value <= 2]*

All other values are only calculated when 2 trace file are analysed in parallel: 

- 'Discrepancy' between 2 chains is calculated as described in the [Phylobayes manual](http://megasun.bch.umontreal.ca/People/lartillot/www/phylobayes3.3e.pdf) and therefore equivalent to the `rel_diff` values that `tracecomp` generates. *[recommended value <= 0.3]*

- 'Gelman & Rubin' stands for Gelman and Rubin's convergence diagnostic and is calculated with the `gelman.diag` function of [coda](https://cran.r-project.org/web/packages/coda/index.html). This tests for convergence between chains by comparing within and between chain variances. If these are roughly equal, then the chains can be considered converged. *[recommended value <= 1.2]*

**As a visual aid, all values that may indicate lack of stationarity or convergence, or insufficient run lengths are highlighted in red.**
