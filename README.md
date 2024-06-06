# bit24_MDFS_workshop

## MDFS installation instruction:

```R
#uncomment if you want a separate version of MDFS installed in local/directory
#in that case, execute two lines below each time before loading MDFS
#.libPaths()-> default_libPath #save it - after installing/loading new MDFS you can set .libPaths to previous value
#.libPaths("local/directory")
if (!require("remotes")) install.packages("remotes")
remotes::install_github("https://github.com/p100mma/mdfs-r")
```
