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

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>

#ifndef PATH_MAX
#	define PATH_MAX _MAX_PATH
#endif

#define PATH_BUF_SIZE PATH_MAX+1

/* Functions for platforms that are missing them. */

/*******************************************************************************
 * Borrowed from rsync
 * Copyright (C) 1998 Andrew Tridgell
 * Copyright (C) 2002 Martin Pool
 * Copyright (C) 2004, 2005, 2006 Wayne Davison
 */

#ifndef HAVE_STRLCPY
size_t strlcpy(char *d, const char *s, size_t bufsize)
{
	size_t len = strlen(s);
	size_t ret = len;
	if (bufsize > 0) {
		if (len >= bufsize)
			len = bufsize-1;
		memcpy(d, s, len);
		d[len] = 0;
	}
	return ret;
}
#endif

#ifndef HAVE_STRLCAT
size_t strlcat(char *d, const char *s, size_t bufsize)
{
	size_t len1 = strlen(d);
	size_t len2 = strlen(s);
	size_t ret = len1 + len2;

	if (len1 < bufsize - 1) {
		if (len2 >= bufsize - len1)
			len2 = bufsize - len1 - 1;
		memcpy(d+len1, s, len2);
		d[len1+len2] = 0;
	}
	return ret;
}
#endif

/*
 * Copyright (c) 1997 Todd C. Miller <Todd.Miller@courtesan.com>
 * Copyright (c) 2009 Reed A. Cartwright
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

char * basename(const char *path)
{
	static char *bname = NULL;
	const char *endp, *startp;

 	if (bname == NULL) {
		bname = (char *)malloc(MAXPATHLEN);
		if (bname == NULL)
			return(NULL);
	}

	/* Empty or NULL string gets treated as "." */
	if (path == NULL || *path == '\0') {
		(void)strcpy(bname, ".");
		return(bname);
	}

	/* Strip trailing slashes */
	endp = path + strlen(path) - 1;
	while (endp > path && *endp == '/')
		endp--;

	/* All slashes becomes "/" */
	if (endp == path && *endp == '/') {
		(void)strcpy(bname, "/");
		return(bname);
	}

	/* Find the start of the base */
	startp = endp;
	while (startp > path && *(startp - 1) != '/')
		startp--;

	if (endp - startp + 2 > MAXPATHLEN) {
		errno = ENAMETOOLONG;
		return(NULL);
	}
	(void)strncpy(bname, startp, endp - startp + 1);
	bname[endp - startp + 1] = '\0';
	return(bname);
}

/******************************************************************************
 * Borrowed from libiberty
 * Public domain
 */

char * basename (const char *name) {
	const char *base;
#if defined (HAVE_DOS_BASED_FILE_SYSTEM)
	/* Skip over the disk name in MSDOS pathnames. */
	if (isalpha(name[0]) && name[1] == ':') 
		name += 2;
#endif
	for (base = name; *name; name++) {
		if (IS_DIR_SEPARATOR (*name)) {
			base = name + 1;
		}
	}
	return (char *) base;
}

