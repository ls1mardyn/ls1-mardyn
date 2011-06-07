#ifndef COMPILE_INFO_H_
#define COMPILE_INFO_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define MAX_INFO_STRING_LENGTH 1024



int get_compiler_info(char **info_str) {
	/* For compiler predfined macros used for identification see
	 * http://predef.sourceforge.net/precomp.html
	 */

	/* Cray compiler */
#if defined(_CRAYC)
	sprintf( *info_str, "Cray %d.%d", _RELEASE, _RELEASE_MINOR);
#endif

	/* GNU compiler */
#if defined(__GNUC__)
# if defined(__GNUC_PATCHLEVEL__)
	sprintf(*info_str, "GNU %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
# else
	sprintf(*info_str, "GNU %d.%d", __GNUC__, __GNUC_MINOR__);
# endif

	/* Intel compiler */
#elif defined(__INTEL_COMPILER)
	int version = (int) __INTEL_COMPILER / 100;
	int revision = ((int) __INTEL_COMPILER % 100) / 10;
	int patch = (int) __INTEL_COMPILER % 10;
	sprintf(*info_str, "Intel %d.%d.%d", version, revision, patch);

	/* PGI compiler */
#elif defined(__PGI)
	sprintf(*info_str, "PGI");

	/* Sun compiler */
#elif defined(__SUNPRO_C)
	int version = (__SUNPRO_C >> 12) & 0xf;
	int revision_digit1 = (__SUNPRO_C >> 8) & 0xf;
	int revision_digit2 = (__SUNPRO_C >> 4) & 0xf;
	int patch = (__SUNPRO_C >> 0 ) & 0xf;
	sprintf(*info_str, "Sun %d.%d%d.%d", version, revision_digit1, revision_digit2, patch);

	/* unknown */
#else
	sprintf(*info_str, "unknown");
#endif

	return 0;
}

#if defined(MPI_VERSION)

int get_mpi_info(char **info_str) {
	/* Intel MPI */
	/* FIXME */

	/* MVAPICH2 */
#if defined(MVAPICH2)
	sprintf(*info_str, "MVAPICH2");

	/* Open MPI */
#elif defined(OPEN_MPI)
	sprintf(*info_str, "Open MPI");

	/* unknown */
#else
	sprintf(*info_str, "unknown");
#endif

	return 0;
}

#endif

int get_compile_time(char **info_str) {
	sprintf(*info_str, "%s %s", __DATE__, __TIME__);
	return 0;
}

int get_timestamp(char **info_str) {
	time_t t;
	struct tm *ts;
	t = time(NULL);
	ts = localtime(&t);
	sprintf(*info_str, "D%d-%d-%d T%d:%d:%d",
	        ts->tm_mday, ts->tm_mon + 1, ts->tm_year + 1900,
	        ts->tm_hour, ts->tm_min, ts->tm_sec);
	return 0;
}

int get_host(char **info_str) {
	char hostname[1024];
	hostname[1023] = '\0';
	gethostname(hostname, 1023);
	sprintf(*info_str, "%s", hostname);
	return 0;
}

std::string getCompileFlags() {
	std::stringstream flags;

#ifdef ENABLE_MPI
	flags << " ENABLE_MPI";
#endif

#ifdef NDEBUG
	flags << " NDEBUG ";
#endif

#ifdef VTK
	flags << " VTK ";
#endif

#ifdef UNIT_TESTS
	flags << " UNIT_TESTS ";
#endif

	return flags.str();
}

#endif /*COMPILE_INFO_H_*/
