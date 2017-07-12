.onLoad=function(libname, pkgname)
{
library.dynam("relateAdmix", pkgname, libname)
}
.onUnload=function(libpath)
{
library.dynam.unload("relateAdmix", libpath)
}