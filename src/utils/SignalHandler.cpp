#include "SignalHandler.h"

#ifndef NDEBUG
#ifdef __GLIBC__
#include <execinfo.h>
#endif
#endif

#include "mardyn_assert.h"

#ifndef NDEBUG
std::string getStackTrace() {
	std::ostringstream ss;
#ifdef __GLIBC__
	size_t size = 10;
	void* array[size];
	size = backtrace(array, size);
	char** symbols = backtrace_symbols(array, size);
	if (symbols != nullptr) {
		ss << "Stack trace:\n";
		for (int i = 0; i < size; ++i)
			ss << "  " << symbols[i] << '\n';
		free(symbols);
	}
#endif
	return ss.str();
}
#endif

std::atomic<int> SignalHandler::_signalBitmaskAtomic = SignalHandler::SIG_NONE;

void SignalHandler::handleSignal(int signalReceived) {
#ifdef ENABLE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	Log::global_log->info()
#ifdef ENABLE_MPI
		<< "[Rank #" << rank << "] "
#endif
		<< "Received signal: " << signalReceived << std::endl;
	int signalBitmask = 0;
	std::ostringstream ossError;
	switch (signalReceived) {
		case SIGINT:
		case SIGTERM:
			signalBitmask = SignalHandler::SIG_STOP;
			break;
		case SIGUSR1:
			signalBitmask = SignalHandler::SIG_USR1;
			break;
#ifndef NDEBUG
		case SIGSEGV:
			ossError << "Segmentation fault (" << signalReceived << ")!\n" << getStackTrace();
			MARDYN_EXIT(ossError.str());
#endif
		default:
			ossError << "Handler caught wrong signal: " << signalReceived;
			MARDYN_EXIT(ossError.str());
	}
	_signalBitmaskAtomic.fetch_or(signalBitmask, std::memory_order_relaxed);
}

SignalHandler* SignalHandler::_enabledInstance = nullptr;

void SignalHandler::enable() {
	if (isEnabled()) {
#ifndef NDEBUG
		Log::global_log->warning() << "Signal handler is already enabled" << std::endl;
#endif
		return;
	}
	if (_enabledInstance != nullptr) {
#ifndef NDEBUG
		Log::global_log->warning() << "Another signal handler was still enabled" << std::endl;
#endif
		_enabledInstance->disable();
	}
	Log::global_log->info() << "Installing signal handler" << std::endl;
	struct sigaction sa {};
	sa.sa_handler = SignalHandler::handleSignal;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = 0;

	sigaction(SIGINT, &sa, &_oldSigInt);
	sigaction(SIGTERM, &sa, &_oldSigTerm);
	sigaction(SIGUSR1, &sa, &_oldSigUsr1);
#ifndef NDEBUG
	sigaction(SIGSEGV, &sa, &_oldSigSegv);
#endif
	_enabledInstance = this;
}

void SignalHandler::disable() {
	if (!isEnabled()) {
#ifndef NDEBUG
		Log::global_log->warning() << "Signal handler is not enabled" << std::endl;
#endif
		return;
	}
	Log::global_log->info() << "Restoring old signal handlers" << std::endl;
	sigaction(SIGINT, &_oldSigInt, nullptr);
	sigaction(SIGTERM, &_oldSigTerm, nullptr);
	sigaction(SIGUSR1, &_oldSigUsr1, nullptr);
#ifndef NDEBUG
	sigaction(SIGSEGV, &_oldSigSegv, nullptr);
#endif
	_enabledInstance = nullptr;
}
