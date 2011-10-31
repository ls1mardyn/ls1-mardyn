#ifndef DOMAIN_DECOMP_TYPES_H
#define DOMAIN_DECOMP_TYPES_H

/** The DomainDecompType datatype provides IDs for the different domain decomposition implementations. */
typedef enum {
	UNKNOWN_DECOMPOSITION = 0, /**< unknown domain decomposition type */
#ifdef ENABLE_MPI
	DOMAIN_DECOMPOSITION, /**< simple domain decomposition */
	KD_DECOMPOSITION, /** KD based domain decomposition */
#endif
	DUMMY_DECOMPOSITION /** Dummy decomposition representing serial program */
} DomainDecompType;

#endif /* DOMAIN_DECOMP_TYPES_H */
