.onAttach <- function(libname, pkgname) {
  
  github_note <- paste0(
    "If you installed Rlibeemd from CRAN, consider installing again from ",
    "R-universe for a version supporting parallelization: ",
    "install.packages('Rlibeemd', repos = 'https://helske.r-universe.dev')")
  
  citation_note <- paste0("\nPlease cite Rlibeemd in publications by using: \n",
    "Luukko PJ, Helske J, R\U00E4s\U00E4nen E (2016). Introducing libeemd:", 
    "A program package for performing the ensemble empirical mode decomposition.",
    "Computational Statistics, 31(2), 545-557. doi: 10.1007/s00180-015-0603-9.")
  
  packageStartupMessage(paste(github_note, citation_note, sep = "\n"))
}