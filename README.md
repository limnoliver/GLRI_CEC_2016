# MMSD_trends

## About this workflow

This project uses [`remake`](https://github.com/richfitz/remake), a file dependency manager that ensures that our analysis scripts get run in a sensible order. from richfitz:

You describe the beginning, intermediate and end points of your analysis, and how they flow together.

* "**targets**" are the points in your analysis.  They can be either files (data files for input; tables, plots, knitr reports for output) or they can be R objects (representing processed data, results, fitted models, etc).
* "**rules**" are how the targets in your analysis relate together and are simply the names of R functions.
* "**dependencies**" are the targets that need to already be made before a particular target can run (for example, a processed data set might depend on downloading a file; a plot might depend on a processed data set).

## Setup

Configuration to use `remake` with this R project is really simple!:

1. install `remake`:
```r
devtools::install_github("richfitz/remake")
```

2. install any missing packages that this project requires that you don't already have:
```r
remake::install_missing_packages()
```


## Building the project

Build this project, or pieces of it, using `remake`. 
```r
library(remake)
# build the entire project:
make() 

# build only one stage of the project:
make(remake_file = "10_load_data.yml")

# build only one target of the project:
make("20_merge_data/doc/progress.csv")
```

## debugging in R

`remake` is "unapologetically R focussed", and supports debugging of functions by inserting `browser()` commands inline to your functions


Want to do things the old fashioned way? you can create the script that `remake` would execute if all targets were out of date:
```r
remake::make_script()
```

If you would like a particular target called outside the main workflow, you can call it directly with `remake::make`:

```r
merged_data <- remake::make("merged_data")
```

## starting fresh
like `make`, you can start a "clean" build:
```r
remake::make("clean")

```
Note that the above command deletes files and also gets rid of R objects. 
Alternatively, you can delete individual targets:
```r
remake::delete("20_merge/doc/progress.csv")
```

### What happens in a build

Subfolders named 'out' and 'log' exist within each numbered folder, and there are a few 'doc' subfolders here and there. On GitHub, these are empty except for README.md files. The README.md files serve as placeholders so that the directories can be versioned and don't need to be created by the project scripts. When you build the project, these folders become populated with data files, figures, etc. ('out'), and ancillary documentation ('doc'). 

## R scripts

What's going on?

### 10_load_data

Raw data files are saved on a private S3 bucket. The function in this step assumes you have a "default" credential set up on your computer. Then, the files are simply downloaded to the "1_get_raw_data/out" folder.



### Dependency tree

[![Dependencies](6_model/remake_diagram.png)

The procedure for making a remake dependency diagram:
```r
remake::diagram()
```
which also takes the same arguments as `remake::make()`, so you can build a diagram for a stage, the whole project, or a single target.

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey (USGS), an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at <https://www.usgs.gov/visual-id/credit_usgs.html#copyright>

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."

[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)](http://creativecommons.org/publicdomain/zero/1.0/)