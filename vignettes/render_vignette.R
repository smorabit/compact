suppressWarnings({
  library(rmarkdown, quietly = TRUE)
  render('simulation_tutorial.Rmd', output_format = 'html_vignette')
})
