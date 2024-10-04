# Splinets 1.5.1

  -  The titles in the dyadic plot of bases are The nth order splinet and they appear at each level of the dyadic structure.  It is better that the user defines the title. There should be one common title over the entire dyadic structure (or, alternatively, different titles indicating which level of dyadic structure one sees). 
The following function in AuxFun.R should be modified:  
`plot.basis`, `plot.obasis`, `plot.splinet`. The above was corrected by completely removing the titles, all `main= ...`. Generally, titles should be avoided in the graphs. The figures captions are the way of describing what is on the graphs. 


  - The scaling of the plot of b-splines in the standard plot (not dyadic one) is not properly scaled.  It is fixed in the version 1.5.1. the above was corrected by changing in `AuxFun` the following function: plot.basis. Before there was `ylim=range(y)`, which included also the first column. It was also a problem in the function plot.obasis() and was fixed as well. 


  - In the manual is: `projedt$sp` should be `poject$sp`. This is fixed in Version 1.5.1


  - Some problems on using 'projection' to get refined representation of a spline:
This is fixed in Version 1.5.1


  -  `print.class()` in `AuxFun.R` file has now added @' export and a more complete description that also appears in the manual. A simple examples are given as well. I could not put it outside of `AuxFun.R` which I would prefer. I could not build the package (maybe cleaning all .Rd files would help?). If one knows how to place it together with all function in the manual, please, do it. This can be a method rather than a function, see `is.splinets()` or 'plot' how this could be done. 


  -  A new argument to `is.splinets()` is given to control the accuracy the conditions in the matrix of the derivatives to be checked. Before it was given by the field in `Splinets`-object, now it can be overiden by 'epsilon' argument of `is.splinets()`. 
To address something that was complained by a reviewer of R-paper. 
