/************************************************************************* 
 * Copyright (c) 2009 Reed A. Cartwright                                 *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *************************************************************************/
#pragma once
#ifndef COMPAT_H
#define COMPAT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stddef.h>

#ifndef HAVE_STRLCPY
size_t strlcpy(char *d, const char *s, size_t bufsize);
#endif
#ifndef HAVE_STRLCAT
size_t strlcat(char *d, const char *s, size_t bufsize);
#endif
#ifndef HAVE_DIRNAME
char * dirname(const char *path);
#endif
#ifndef HAVE_BASENAME
char * basename(const char *path);
#endif

#ifndef DIR_SEPARATOR
#define DIR_SEPARATOR '/'
#endif

#if defined (_WIN32) || defined (__MSDOS__) || defined (__DJGPP__) || \
  defined (__OS2__)
#define HAVE_DOS_BASED_FILE_SYSTEM
#ifndef DIR_SEPARATOR_2 
#define DIR_SEPARATOR_2 '\\'
#endif
#endif

/* Define IS_DIR_SEPARATOR.  */
#ifndef DIR_SEPARATOR_2
# define IS_DIR_SEPARATOR(ch) ((ch) == DIR_SEPARATOR)
#else /* DIR_SEPARATOR_2 */
# define IS_DIR_SEPARATOR(ch) \
	(((ch) == DIR_SEPARATOR) || ((ch) == DIR_SEPARATOR_2))
#endif /* DIR_SEPARATOR_2 */

#ifndef DIR_SEPARATOR_STR
#	ifdef HAVE_DOS_BASED_FILE_SYSTEM
#		define DIR_SEPARATOR_STR "\\"
#	else
#		define DIR_SEPARATOR_STR "/"
#	endif
#endif

#endif
