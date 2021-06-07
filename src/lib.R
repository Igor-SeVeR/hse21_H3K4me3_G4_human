# At first we set working directory

setwd('C:/minorBio/final_project/src')

# libraries

# at first, we install all libraries, which are needed to make the project

# method to check, if library is already installed and if not -> install
package.check <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
}

package.check("ggplot2")
package.check("dplyr")
package.check("tidyr")
package.check("tibble")

# constants

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'