.onLoad <- function(libname, pkgname){
  load(system.file("data/vsearch_executable.rds", package = "Rsearch"))
  options(Rsearch.vsearch_executable = vsearch_executable)
}
