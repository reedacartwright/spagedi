#pragma once
#ifndef SPAGEDI_CONFIG_H
#define SPAGEDI_CONFIG_H

#define SPAGEDI_VERSION "@SPAGEDI_VERSION@"

#cmakedefine  HAVE_UNISTD_H 1
#cmakedefine  HAVE_DIRECT_H 1
#cmakedefine  HAVE_LIBGEN_H 1
#cmakedefine  HAVE_CONIO_H 1

#cmakedefine HAVE_CHDIR 1
#cmakedefine HAVE_STRLCPY 1
#cmakedefine HAVE_STRLCAT 1
#cmakedefine HAVE_BASENAME 1
#cmakedefine HAVE_DIRNAME 1
#cmakedefine HAVE__SPLITPATH 1

#endif
