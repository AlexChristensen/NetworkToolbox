.onload <- function(libname, pkgname)
{library.dynam("NetworkToolbox",package=pkgname,lib.loc=libname)}

.onAttach <- function(libname, pkgname)
{
    temp <- packageDescription("NetworkToolbox")
    msg <- paste("Package: ",temp$Package,": ",temp$Title,"\n",
               "Version: ",temp$Version,"\n",
               "Updated on: ",
               temp$Date,"\n", sep="")
    msg <- paste(msg,"Maintainer: Alexander P. Christensen, University of North Carolina at Greensboro\n",sep="")
    msg <- paste(msg,"Contributors: Guido Previde Massara, University College London\n",sep="")
    msg <- paste(msg,'For citation information, type citation("NetworkToolbox")\n')
    msg <- paste(msg,'For vignettes, see <https://doi.org/10.32614/RJ-2018-065> \n')
    msg <- paste(msg,"For bugs and errors, submit an issue to <https://github.com/AlexChristensen/NetworkToolbox/issues>")
    packageStartupMessage(msg)
}