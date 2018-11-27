#ifndef PACKAGE_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#ifdef PACKAGE_STRING
#	define VERSION PACKAGE_STRING
#else
#	define VERSION "SPAGeDi " SPAGEDI_VERSION
#endif

#endif

