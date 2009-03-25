#cmakedefine  CURSES_HAVE_CURSES_H 1
#cmakedefine  CURSES_HAVE_NCURSES_H 1
#cmakedefine  CURSES_HAVE_NCURSES_NCURSES_H 1
#cmakedefine  CURSES_HAVE_NCURSES_CURSES_H 1
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

/******************************************************************************
 * Borrowed from libiberty
 * Public domain
 */
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

