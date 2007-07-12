#include "utils/Log.h"

utils::CommandLineLogger utils::Log::_logger;


utils::Log::Log(const std::string& className):
  _className( className ) {
}

utils::Log::~Log() {
}

void utils::Log::debug(const std::string& methodName, const std::string& message) const {
  Log::_logger.debug( _className, methodName, message);
}

void utils::Log::info(const std::string& methodName, const std::string& message) const {
  _logger.info( _className, methodName, message);
}

void utils::Log::warning(const std::string& methodName, const std::string& message) const {
  _logger.warning( _className, methodName, message);
}

void utils::Log::error(const std::string& methodName, const std::string& message) const {
  _logger.error( _className, methodName, message);
}
