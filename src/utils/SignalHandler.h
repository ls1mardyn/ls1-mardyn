#ifndef SIGNAL_HANDLER_H_
#define SIGNAL_HANDLER_H_

#include <atomic>
#include <csignal>

/**
 * @brief Custom handler for SIGINT, SIGTERM, SIGUSR1, (and SIGSEGV for debug builds)
 */
class SignalHandler {
public:
	/**
	 * @brief Replace the current signal handlers with this signal handler
	 */
	void enable();

	/**
	 * @brief Check if this signal handler is currently being used
	 */
	inline bool isEnabled() { return _enabledInstance == this; };

	/**
	 * @brief Disable the signal handler and restore old signals
	 */
	void disable();

	/**
	 * @brief Update the signals bitmask with received signals across all ranks from the last interval since calling
	 * this function.
	 */
	inline void syncReceivedSignals() {
		_signalBitmask = _signalBitmaskAtomic.exchange(SIG_NONE, std::memory_order_relaxed);
#ifdef ENABLE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &_signalBitmask, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);
#endif
	}

	/**
	 * @brief Return the bitmask of received signals
	 * @returns A bitmask of \ref SignalBits
	 */
	inline int getBitmask() { return _signalBitmask; };

	/**
	 * @brief Bitmask mapping for individual signals that can be handled
	 */
	enum SignalBits {
		SIG_NONE = 0,		// No signal
		SIG_STOP = 1 << 0,	// SIGINT/SIGTERM
		SIG_USR1 = 1 << 1,	// SIGUSR1
	};

	SignalHandler() {};
	~SignalHandler() {
		if (isEnabled())
			disable();
	};

	// Prevent accidental creation of multiple signal handlers
	SignalHandler(const SignalHandler&) = delete;
	SignalHandler& operator=(const SignalHandler&) = delete;
	SignalHandler(SignalHandler&&) = delete;
	SignalHandler& operator=(SignalHandler&&) = delete;

private:
	// Not a singleton, but avoiding undefined behaviour throuh multiple signal handlers
	static SignalHandler* _enabledInstance;
	// Store old signal handlers and restore them later
	struct sigaction _oldSigInt;
	struct sigaction _oldSigTerm;
	struct sigaction _oldSigUsr1;
	struct sigaction _oldSigSegv;
	/**
	 * @brief Bitmask of received signals (Changed asynchronously)
	 */
	static std::atomic<int> _signalBitmaskAtomic;
	/**
	 * @brief Stores signals received on any rank (Handling on all ranks is deferred to \ref syncReceivedSignals)
	 */
	static void handleSignal(int signalReceived);
	/**
	 * @brief Copy of received signals (Consistent across ranks)
	 */
	int _signalBitmask = SIG_NONE;
};

#endif	// SIGNAL_HANDLER_H_
