/*
 * SigsegvHandler.h
 *
 *  Created on: May 11, 2017
 *      Author: seckler
 */
#pragma once

#ifdef ENABLE_SIGHANDLER

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include "mardyn_assert.h"
#include "utils/mardyn_assert.h"


void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  MARDYN_EXIT(1);
}

void registerSigsegvHandler() {
  signal(SIGSEGV, handler);   // install our handler
}

#endif

