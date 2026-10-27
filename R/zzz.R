.onAttach <- function(libname, pkgname)
{
    temp <- packageDescription("NetworkToolbox")
    msg <- paste("Package: ",temp$Package,": ",temp$Title,"\n",
               "Version: ",temp$Version,"\n",
               "Updated on: ",
               temp$Date,"\n", sep="")
    msg <- paste(msg,'For citation information, type citation("NetworkToolbox")\n')
    msg <- paste(msg,'For vignettes, see <https://doi.org/10.32614/RJ-2018-065> \n')
    msg <- paste(msg,"For bugs and errors, submit an issue to <https://github.com/AlexChristensen/NetworkToolbox/issues>")
    packageStartupMessage(msg)
}