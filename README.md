# bit24_MDFS_workshop

Software needed:

- `R >= 3.4.0`
- `Rtools` - **we are going to use a development version of MDFS package** which has to be installed from source. Newest version of MDFS currently available on CRAN does **not** contain all the functionalities we wish to present.
https://cran.r-project.org/bin/windows/Rtools

## MDFS installation instruction:

```R
#uncomment if you want a separate version of MDFS installed in local/directory
#in that case, execute two lines below each time before loading MDFS
#.libPaths()-> default_libPath #save it - after installing/loading new MDFS you can set .libPaths to previous value
#.libPaths("local/directory")
if (!require("remotes")) install.packages("remotes")
remotes::install_github("https://github.com/p100mma/mdfs-r")
```
