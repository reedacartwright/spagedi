#ifndef PACKAGE_H

#ifdef HAVE_VERSION_H
#	include "version.h"
#endif

#ifdef PACKAGE_STRING
#	define VERSION PACKAGE_STRING
#else
#	define  VERSION "SPAGeDi 1.4a (build 11-01-2013)"
#endif

#endif

