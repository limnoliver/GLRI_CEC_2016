render_html <- function(filename_md, output_file) {
  
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  rmarkdown::render(filename_md, 
                    output_file = basename(output_file),
                    output_dir = dirname(output_file),
                    output_format = "html_document", 
                    quiet=TRUE)
}