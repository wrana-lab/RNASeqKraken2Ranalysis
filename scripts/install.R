renv::init(bare = TRUE)
renv::clean()
renv::status()

install.packages("here", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("tidyr", dependencies = TRUE)
install.packages("stringr", dependencies = TRUE)
install.packages("purrr", dependencies = TRUE)
install.packages("tibble", dependencies = TRUE)
install.packages("ggrepel", dependencies = TRUE)
install.packages("kableExtra", dependencies = TRUE)
install.packages("gridExtra", dependencies = TRUE)

install.packages("foreach")
install.packages("doParallel")

renv::snapshot()
renv::restore()
renv::diagnostics(project = NULL)

# Temporary
install.packages("argparse")
install.packages("reshape2")
install.packages("rsthemes")
install.packages("viridis")

packageVersion("here")
packageVersion("ggplot2")
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("purrr")




