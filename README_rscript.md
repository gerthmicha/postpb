#### pbpplot.R

This is the standalone version to be run as executable R script

```
Usage: /usr/bin/pbplot.R [options]

Options:
	-b NUMBER, --burnin=NUMBER
		number of burnin iterations to discard [required]

	-f, --file
		Shall plot be saved to pdf file instead of disp? [default: FALSE]

	-h, --help
		Show this help message and exit

NOTE:	This script plots the parameters estimated during a phylobayes run.
	It should be run in a directory containing 1–4 trace files created by phylobayes.
	Requires R packages ggplot2, optparse, and gridExtra.
```
