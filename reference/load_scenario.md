# Load a scenario folder

Reads the input workbooks (combined data, optional dummy data, model
parameter sheet) from a scenario directory and assembles the inputs
needed by \[fitCompartmentalModel()\] / \[build_compartmental_model()\],
plus an output plot path.

## Usage

``` r
load_scenario(
  scenario_dir,
  combined_file = "dataCombined.xlsx",
  dummy_file = "dataDummy.xlsx",
  params_file = "modelParams.xlsx",
  plots_subdir = "Plots",
  plot_file = "Plot.pdf"
)
```

## Arguments

- scenario_dir:

  Folder containing the scenario inputs.

- combined_file:

  Combined-data workbook name.

- dummy_file:

  Optional dummy-data workbook name.

- params_file:

  Model parameter sheet name.

- plots_subdir:

  Sub-folder for output plots.

- plot_file:

  Default output plot file name.

## Value

A list with the loaded \`dataCombined\`, \`dataDummy\`, \`modelParams\`,
and the resolved plot path.

## Examples

``` r
# The package ships a minimal example scenario as csv:
mini <- system.file("extdata", "minimal", package = "compfit")
sc <- load_scenario(mini,
                    combined_file = "dataCombined.csv",
                    dummy_file    = "dataDummy.csv",
                    params_file   = "modelParams.csv")
str(sc, max.level = 1)
#> List of 6
#>  $ dir         : chr "/home/runner/work/_temp/Library/compfit/extdata/minimal"
#>  $ dataCombined: spc_tbl_ [2 × 11] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
#>   ..- attr(*, "spec")=
#>   .. .. cols(
#>   .. ..   .default = col_character(),
#>   .. ..   Label = col_character(),
#>   .. ..   Formula = col_character(),
#>   .. ..   Likelihood = col_character(),
#>   .. ..   Weight = col_character(),
#>   .. ..   Average = col_character(),
#>   .. ..   `2015` = col_character(),
#>   .. ..   `2016` = col_character(),
#>   .. ..   `2017` = col_character(),
#>   .. ..   `2018` = col_character(),
#>   .. ..   `2019` = col_character(),
#>   .. ..   `2020` = col_character()
#>   .. .. )
#>   ..- attr(*, "problems")=<externalptr> 
#>  $ dataDummy   : spc_tbl_ [1 × 2] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
#>   ..- attr(*, "spec")=
#>   .. .. cols(
#>   .. ..   Label = col_character(),
#>   .. ..   Formula = col_character()
#>   .. .. )
#>   ..- attr(*, "problems")=<externalptr> 
#>  $ modelParams : spc_tbl_ [4 × 10] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
#>   ..- attr(*, "spec")=
#>   .. .. cols(
#>   .. ..   .default = col_character(),
#>   .. ..   Level_1 = col_character(),
#>   .. ..   Others = col_character(),
#>   .. ..   States = col_character(),
#>   .. ..   Functions = col_character(),
#>   .. ..   Parameters = col_character(),
#>   .. ..   Conditions = col_character(),
#>   .. ..   Linear1 = col_character(),
#>   .. ..   Quadratic1 = col_character(),
#>   .. ..   Linear2 = col_character(),
#>   .. ..   Quadratic2 = col_character()
#>   .. .. )
#>   ..- attr(*, "problems")=<externalptr> 
#>  $ plot_dir    : chr "/home/runner/work/_temp/Library/compfit/extdata/minimal/Plots"
#>  $ plot_path   : chr "/home/runner/work/_temp/Library/compfit/extdata/minimal/Plots/Plot.pdf"
```
