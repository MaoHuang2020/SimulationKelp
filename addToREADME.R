#### writeREADME
# This is a function to write lines to a README.md file.  In .md, lines have to end with two spaces to cause a carriage return.
addToREADME <- function(strVec, append=T){
  # Starting fresh: point user to correct README.md
  if (!append) write_lines("[See README.md in docs folder](./README.md)", "README.md", sep="  \n", append=F)
  # Add text with two space to get carriage return
  write_lines(c(strVec, ""), "./README.md", sep="  \n", append=append)
}
