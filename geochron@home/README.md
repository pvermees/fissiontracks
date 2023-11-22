### Summary

**[geochron@home](https://github.com/pvermees/geochron-at-home)** is a
web platform that serves as both an archiving tool for published
fission track data, and a crowd-sourcing platform for new fission
track images.

**FT2GaH.R** is an **R** script to convert TIF-images and XML-files
produced by AutoScan's FastTrack software to JPG-stacks and JSON-files
that are suitable for uploading to **geochron@home**.

### prerequisites

1. **R**, from [CRAN](https://r-project.org)

2. The **XML** and **jsonlite** packages. At the **R** console:

```
install.packages('XML')
install.packages('jsonlite')
```

3. **ImageMagick**, from [https://imagemagick.org](https://imagemagick.org)

4. **FT2GaH.R**, from this GitHub repository: [FT2GaH.R](https://github.com/pvermees/fissiontracks/tree/master/geochron@home/FT2GaH.R)

### Example

```
source("/path/to/FT2GaH.R")
FT2GaH(idir="/path/to/FT/data",odir="/path/to/GaH/data",
       xml="results.xml",sample="samplename")
```

where

1. `idir` expects the full path to your FastTrack data input folder;

2. `odir` expects the full path to your geochron@home output folder;

3. `xml` expects the name of the XML results file (which should be
located inside the input folder;

4. `sample` expects the sample name.

On Windows computers, you may need to provide the full path to
ImageMagick's `convert` function via the optional `IMpath`
argument. For example:

```
source("/path/to/FT2GaH.R")
FT2GaH(idir="/path/to/FT/data",odir="/path/to/GaH/data",
       xml="results.xml",sample="samplename",
       IMpath="\"C:\\Program Files\\ImageMagick\"")
```

Note the necessary `\` escape characters.