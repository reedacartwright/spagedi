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

#include "compat.h"

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#ifndef PATH_MAX
#	define PATH_MAX _MAX_PATH
#endif

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

/******************************************************************************/

/* Treat the terminal null as a sep as well to simplify loop. */
#define IS_DIR_SEPARATOR_2(x) ((IS_DIR_SEPARATOR(x)) || ((x) == '\0'))

#ifndef HAVE_BASENAME
char * basename(const char *path) {
	static char bname[PATH_MAX];
	const char *p, *endp, *startp;

	/* Empty or NULL string gets treated as "." */
	if (path == NULL || *path == '\0') {
		strcpy(bname, ".");
		return bname;
	}
#if defined (HAVE_DOS_BASED_FILE_SYSTEM)
	/* Skip over the disk name in MSDOS pathnames. */
	if (isalpha(path[0]) && path[1] == ':') {
		path += 2;
	/* Check for NULL string */
		if (*path == '\0') {
			strcpy(bname, ".");
			return bname;
		}
	}
#endif
	startp = p = path;
	endp = path+1;
	for(p=path;*p;++p) {
		if(IS_DIR_SEPARATOR(*p)) {
			if(!IS_DIR_SEPARATOR_2(*(p+1)))
				/* SEP NSEP = start of new dir */
				startp = p+1;
		} else if(IS_DIR_SEPARATOR_2(*(p+1)))
			/* NSEP SEP = end of new dir */
			endp = p+1;
	}
	if(endp-startp+1 > PATH_MAX)
		return NULL;
	strncpy(bname, startp, endp-startp);
	bname[endp-startp] = '\0';
	return bname;
}
#endif

#ifndef HAVE_DIRNAME
char * dirname(const char *path) {
	static char dname[PATH_MAX];
	const char *p, *startp, *endp, *endpp;

	/* Empty or NULL string gets treated as "." */
	if (path == NULL || *path == '\0') {
		strcpy(dname, ".");
		return dname;
	}
	startp = path;
#if defined (HAVE_DOS_BASED_FILE_SYSTEM)
	/* Skip over the disk name in MSDOS pathnames. */
	if (isalpha(startp[0]) && startp[1] == ':') {
		startp += 2;
	/* Check for NULL string */
		if (*startp == '\0') {
			/* if PATH_MAX > 3 this will work */
			strcpy(dname, path);
			strcat(dname, ".");
			return dname;
		}
	}
#endif
	endp = endpp = p = startp;
	for(;*p;++p) {
		if(IS_DIR_SEPARATOR(*p)) {
			if(!IS_DIR_SEPARATOR_2(*(p+1)))
				/* SEP NSEP = start of new dir */
				endp = endpp;
		} else if(IS_DIR_SEPARATOR_2(*(p+1)))
			/* NSEP SEP = end of new dir */
			endpp = p+1;
	}
	/* no dir? */
	if(endp == startp) {
		p = IS_DIR_SEPARATOR(*endp) ? DIR_SEPARATOR_STR : ".";
		/* check to see if we need to copy a drive letter */
#ifdef HAVE_DOS_BASED_FILE_SYSTEM
		if(startp > path) {
			strncpy(dname, path, startp-path);
			dname[startp-path] = '\0';
			strcat(dname, p);
		} else
#endif
		strcpy(dname, p);
		return dname;
	}
	if(endp-path+1 > PATH_MAX)
		return NULL;
	strncpy(dname, path, endp-path);
	dname[endp-path] = '\0';
	return dname;
}
#endif
