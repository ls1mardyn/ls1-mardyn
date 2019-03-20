#ifndef COMPILE_INFO_H_
#define COMPILE_INFO_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define MAX_INFO_STRING_LENGTH 1024



int get_compiler_info(char *info_str) {
	/* For compiler predfined macros used for identification see
	 * http://predef.sourceforge.net/precomp.html
	 */

	/* Cray compiler */
#if defined(_CRAYC)
	sprintf( info_str, "Cray %d.%d", _RELEASE, _RELEASE_MINOR);

	/* Intel compiler */
#elif defined(__INTEL_COMPILER)
	int version = (int) __INTEL_COMPILER / 100;
	int revision = ((int) __INTEL_COMPILER % 100) / 10;
	int patch = (int) __INTEL_COMPILER % 10;
	sprintf(info_str, "Intel %d.%d.%d", version, revision, patch);

	/* PGI compiler */
#elif defined(__PGI)
	int version = (int) __PGIC__;
	int revision = (int) __PGIC_MINOR__;
	int patch = (int) __PGIC_PATCHLEVEL__;
	sprintf(info_str, "PGI %d.%d.%d", version, revision, patch);

	/* Sun compiler */
#elif defined(__SUNPRO_C)
	int version = (__SUNPRO_C >> 12) & 0xf;
	int revision_digit1 = (__SUNPRO_C >> 8) & 0xf;
	int revision_digit2 = (__SUNPRO_C >> 4) & 0xf;
	int patch = (__SUNPRO_C >> 0 ) & 0xf;
	sprintf(info_str, "Sun %d.%d%d.%d", version, revision_digit1, revision_digit2, patch);

#elif defined(_SX)
	sprintf(info_str, "NEC SX, rev. %d", __SXCXX_REVISION);

	/* Clang compiler */
#elif defined(__clang__)
	int version = __clang_major__;
	int revision = __clang_minor__;
	int patch = __clang_patchlevel__;
	sprintf(info_str, "Clang %d.%d.%d", version, revision, patch);

	/* GNU compiler */
#elif defined(__GNUC__)
# if defined(__GNUC_PATCHLEVEL__)
	sprintf(info_str, "GNU %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
# else
	sprintf(info_str, "GNU %d.%d", __GNUC__, __GNUC_MINOR__);
# endif
	
	/* unknown */
#else
	sprintf(info_str, "unknown");
#endif

	return 0;
}

int get_mpi_info(char *info_str) {

#ifdef ENABLE_MPI

	#if defined(MPI_VERSION) and defined(MPI_SUBVERSION)

		#if defined(MVAPICH2)
			sprintf(info_str, "with MVAPICH2 (MPI %d.%d)", MPI_VERSION, MPI_SUBVERSION);
		#elif defined(OPEN_MPI)
			sprintf(info_str, "with Open MPI %d.%d.%d (MPI %d.%d)", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION, MPI_VERSION, MPI_SUBVERSION);
		#elif defined(CRAY_MPICH_VERSION)
			sprintf(info_str, "with Cray MPI (MPI %d.%d)", MPI_VERSION, MPI_SUBVERSION);
		#elif defined(I_MPI_VERSION)
			sprintf(info_str, "with Intel MPI %s (MPI %d.%d)", I_MPI_VERSION, MPI_VERSION, MPI_SUBVERSION);
		#else
			// unknown
			sprintf(info_str, "with unknown MPI (version %d.%d)", MPI_VERSION, MPI_SUBVERSION);
		#endif

	#else /* MPI_VERSION and MPI_SUBVERSION */
			// I guess this can never happen, as it will fail on #include "mpi.h"?
	#endif /* MPI_VERSION and MPI_SUBVERSION */

#else /* ENABLE_MPI */
	sprintf(info_str, "%s", "without MPI support.");
#endif /* ENABLE_MPI */

	return 0;
}


int get_compile_time(char *info_str) {
	sprintf(info_str, "%s %s", __DATE__, __TIME__);
	return 0;
}

int get_timestamp(char *info_str) {
	time_t t;
	struct tm *ts;
	t = time(NULL);
	ts = localtime(&t);
	sprintf(info_str, "D%d-%d-%d T%d:%d:%d",
	        ts->tm_mday, ts->tm_mon + 1, ts->tm_year + 1900,
	        ts->tm_hour, ts->tm_min, ts->tm_sec);
	return 0;
}

int get_host(char *info_str) {
	char hostname[1024];
#ifdef _SX
	strcpy(hostname, "unknown");
#else
	hostname[1023] = '\0';
	gethostname(hostname, 1023);
#endif
	sprintf(info_str, "%s", hostname);
	return 0;
}

void get_precision_info(char *info_str) {
#if defined(MARDYN_SPSP)
	sprintf(info_str, "%s", "Single (SPSP)");
#elif defined(MARDYN_SPDP)
	sprintf(info_str, "%s", "Mixed (SPDP)");
#else
	sprintf(info_str, "%s", "Double (DPDP)");
#endif
}

void get_intrinsics_info(char *info_str) {
#if VCP_VEC_TYPE==VCP_NOVEC
	sprintf(info_str, "%s", "none");
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	sprintf(info_str, "%s", "SSE3");
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	sprintf(info_str, "%s", "AVX");
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	sprintf(info_str, "%s", "AVX2");
#elif VCP_VEC_TYPE==VCP_VEC_KNL
	sprintf(info_str, "%s", "KNL masking");
#elif VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
	sprintf(info_str, "%s", "KNL gather/scatter");
#elif VCP_VEC_TYPE==VCP_VEC_AVX512F
	sprintf(info_str, "%s", "SKX masking");
#elif VCP_VEC_TYPE==VCP_VEC_AVX512F_GATHER
	sprintf(info_str, "%s", "SKX gather/scatter");
#endif
}

void get_rmm_normal_info(char *info_str) {
#if defined(ENABLE_REDUCED_MEMORY_MODE)
	sprintf(info_str, "%s", "reduced memory mode (RMM). Not all features work in this mode.");
#else
	sprintf(info_str, "%s", "normal mode. All features work in this mode.");
#endif
}

void get_openmp_info(char *info_str) {
#if defined(_OPENMP)
	sprintf(info_str, "%s%d%s", "with OpenMP support, dated ", _OPENMP, ".");
#else
	sprintf(info_str, "%s", "without OpenMP support.");
#endif

}

#endif /*COMPILE_INFO_H_*/
