#!/bin/bash

# check if patch is already applied
VERSION=$(sed --quiet -e 's/#define * VERSION //gp' config.h)

if [[ ${VERSION} =~ .*MARDYN ]]; then
    echo "Quicksched is already patched!"
    exit 0
else
    sed -i -e "s/#define VERSION ${VERSION}/#define QSCHED_VERSION ${VERSION/%\"/_MARDYN\"}/g" config.h
fi

# quicksched has an on function 'error' to print error messages. Since it is not prefixed it interferes with other libraries' error message functions.
sed -i -e 's/error(/qsched_error(/g' src/*

# same with the 'message function'
sed -i -e 's/message(/qsched_message(/g' src/*

# the body of qsched_pthread_run() is missing an if def
sed -i -e 's/void \*qsched_pthread_run ( void \*in ) {/void *qsched_pthread_run ( void *in ) {\n#if defined( HAVE_PTHREAD  )/g' src/qsched.c
# corresponding endif
sed -i -e 'N;s/main loop\. \*\//main loop. *\/\n#endif/g' src/qsched.c

# empty definition of macro TIMER_TOC is missing its arguments
sed -i -e 'N;s/define TIMER_TOC\n/define TIMER_TOC(s,tid)\n/g' src/qsched.h

# fixes undefined reference to getticks
sed -i -e 's/INLINE ticks getticks(void)/static INLINE ticks getticks(void)/g' src/cycle.h

# At the moment only the OpenMP version of Quicksched is supported
sed -i -e 's/#define HAVE_PTHREAD /\/\/#define HAVE_PTHREAD /g' config.h

# sane protection from crashing the compiler by including quicksched twice
sed -i '1i#ifndef QUICKSCHED_H_' src/quicksched.h
sed -i '2i#define QUICKSCHED_H_' src/quicksched.h
sed -i '$a#endif /* QUICKSCHED_H_ */' src/quicksched.h
