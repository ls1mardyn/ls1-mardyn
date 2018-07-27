#!/bin/bash

# check if patch is already applied
VERSION=$(sed --quiet --expression 's/#define * VERSION //gp' config.h)
if [[ -z ${VERSION} ]]; then
    VERSION=$(sed --quiet --expression 's/#define *QSCHED_VERSION //gp' config.h)

    if [[ -z ${VERSION} ]]; then
        echo "Version could not be determined, aborting!"
        exit 1
    fi
fi

if [[ ${VERSION} =~ .*MARDYN ]]; then
    echo "Quicksched is already patched!"
    exit 0
else
    sed --in-place --expression "s/#define VERSION ${VERSION}/#define QSCHED_VERSION ${VERSION/%\"/_MARDYN\"}/g" config.h
fi

# quicksched has an on function 'error' to print error messages. Since it is not prefixed it interferes with other libraries' error message functions.
sed --in-place --expression 's/error(/qsched_error(/g' src/*

# same with the 'message function'
sed --in-place --expression 's/message(/qsched_message(/g' src/*

# the body of qsched_pthread_run() is missing an if def
sed --in-place --expression 's/void \*qsched_pthread_run ( void \*in ) {/void *qsched_pthread_run ( void *in ) {\n#if defined( HAVE_PTHREAD  )/g' src/qsched.c
# corresponding endif
sed --in-place --expression 'N;s/main loop\. \*\//main loop. *\/\n#endif/g' src/qsched.c

# empty definition of macro TIMER_TOC is missing its arguments
sed --in-place --expression 'N;s/define TIMER_TOC\n/define TIMER_TOC(s,tid)\n/g' src/qsched.h

# fixes undefined reference to getticks
sed --in-place --expression 's/INLINE ticks getticks(void)/static INLINE ticks getticks(void)/g' src/cycle.h

# At the moment only the OpenMP version of Quicksched is supported
sed --in-place --expression 's/#define HAVE_PTHREAD /\/\/#define HAVE_PTHREAD /g' config.h

# sane protection from crashing the compiler by including quicksched twice
# also wrapping everything so no ifdefs are necessary in the rest of the code
# this goes directly after the header
sed --in-place --expression 'N;s|\*\*/|**/\
#ifndef QUICKSCHED_H_\
#define QUICKSCHED_H_\
#ifdef QUICKSCHED\
#ifdef __cplusplus\
extern "C"\
#endif\
{|' \
src/quicksched.h
# this goes to the bottom
sed --in-place '$a}\
#endif /* QUICKSCHED */\
#endif /* QUICKSCHED_H_ */' \
src/quicksched.h
