# Close the null PDF device opened in setup.R
if (!interactive() && names(dev.cur()) == "pdf") {
  dev.off()
}
