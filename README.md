# MSsary
*Mass spectrometry in R made easy*
* * *

![Animated ion plot](vignettes/pics/animIon.gif)

## Introduction
Mssary is a package for R under development, with the aim to lessen the burden for both users and developers that needs to interact with raw mass spectrometry data. The scope overlaps somewhat with the popular package xcms, but generally takes a different path in terms of API and infrastructure that makes it easy to work with arbitrarily large datasets while maintaining a quick and efficient link to the underlying raw data.

At the onset MSsary will provide the same algorithms as xcms for deriving grouped peak areas from raw data, but it will provide a rich API for extending the package with additional algorithms. The long term goal is to port all relevant open source algorithms to MSsary (e.g. the OpenMS suite), while also spur the development of new ones directly for MSsary.

MSsary is currently in alpha stage of development, and a stable version from this repository cannot be guarantied. If you have any queries about the project feel free to contact me.
