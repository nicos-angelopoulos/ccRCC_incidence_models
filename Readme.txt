This reproduces the age related results from our paper:
   Timing the Landmark Events in the Evolution of Clear Cell Renal Cell Cancer: TRACERx Renal
   https://doi.org/10.1016/j.cell.2018.02.020

You need to download and unpack the accompanying data in order to run this analysis.

The data should be installed in relative directory data/ with the following sub-directories:

  *  patient_data/  sample_data/  sampler_data/  summary/

To reproduce the knitr generated html file that containes all the code and plots, do:

% 
  R --no-init-file --no-restore-data
> 
   library("knitr");  knit("Chronology_supp_material.Rmd")
>
  library("markdown"); markdownToHTML("Chronology_supp_material.md","Chronology_supp_material.html")

Authors
   Peter Campbell
   Tom Mitchell
   Nicos Angelopoulos
