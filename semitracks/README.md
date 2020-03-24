### Summary

**semitracks.R** is an **R** script to convert semi-track length
  distributions to confined fission track length distributions and
  vice versa. The script is provided with example datasets from the
  University of Melbourne. All the code in this README file can be run
  by simply downloading the contents of
  [http://github.com/pvermees/fissiontracks/semitracks/](http://github.com/pvermees/fissiontracks/semitracks/)
  onto your computer. No other packages are required.

### Examples

1. Forward modelling the anisotropic length distribution of confined
fission tracks and semi-tracks for a three component mixture
containing 20\% tracks of 7 $\pm$ 0.5 $\mu$m long, 40\% tracks of 10
$\pm$ 0.5 $\mu$m long, and 40\% tracks of 14 $\pm$ 0.5 $\mu$m long:


```
source('semitracks.R') # load semitracks.R
P <- c(0.2,0.4,0.4)    # proportions (add up to 1)
M <- c(7,10,14)        # peak locations (microns)
S <- 0.5               # peak width (microns)
plotModel(P,M,S)       # plot the forward model
```

2. Estimating the confined fission track length distribution from a
dataset of 1004 semi-track length measurements:

```
dat <- read.data(fname='track-lengths/UM10-10 Brogo 3D Semi Tracks with Dper.csv',
                 confined=FALSE,skip=5,cols=c(5,9))
fit <- invert(dat,confined=FALSE,ncomp=2)
plotModel(P=fit$P,M=fit$M,S=fit$S,d=dat)
```

where **skip** indicates the length of the header; **cols** marks the
columns containing the length and angle measurements; and **ncomp** is
the number of components.  **fit** is a list with the best fit
proportions (**fit\$P**), modal lengths (**fit\$M**) and peak width
(**fit\$S**).

3. Estimating the semi-track track length distribution from a dataset
of 154 confined track length measurements:

```
dat <- read.data(fname='track-lengths/UM10-10 Brogo 3D Confined Lengths.csv',
                 confined=TRUE,skip=5,cols=c(5,9))
fit <- invert(dat,confined=TRUE,ncomp=2)
plotModel(P=fit$P,M=fit$M,S=fit$S,d=dat)
```

# Further information

Galbraith, R.F., 2005, *Statistics for fission track analysis*. CRC
Press.

Laslett, G.M. and Galbraith, R.F., 1996. Statistical properties of
semi-tracks in fission track analysis. *Radiation measurements*,
26(4), pp.565-576.

Li, Q., Gleadow, A., Seiler, C., Kohn, B., Vermeesch, P., Carter,
A. and Hurford, A., 2018. Observations on three-dimensional
measurement of confined fission track lengths in apatite using digital
imagery. *American Mineralogist*, 103(3), pp.430-440.

# Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

# License

This software is released under the GPL-3 License