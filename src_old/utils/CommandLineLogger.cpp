#include "utils/CommandLineLogger.h"
#include <time.h>
#include <sys/utsname.h>
#include "utils/Globals.h"
#include "utils/LogFilter.h"

#include <cassert>

#ifdef Parallel
#include "parallel/Node.h"
#endif

utils::CommandLineLogger::CommandLineLogger( std::ostream& out ):
  _out( out ) {
  _out.setf( std::ios_base::scientific, std::ios_base::floatfield );
  _out.precision(20);
}

utils::CommandLineLogger& utils::CommandLineLogger::operator=(const CommandLineLogger& rhs) {
  return *this;
}

utils::CommandLineLogger::CommandLineLogger(const CommandLineLogger& param): _out(param._out) {}

utils::CommandLineLogger::~CommandLineLogger() {
  _out.flush();
}

void utils::CommandLineLogger::debug(const std::string& className, const std::string& methodName, const std::string& message) {
  if ( LOG_DEBUG_INFO && WRITE_LOG_INFORMATION ) {
  	
  if ( !utils::LogFilter::writeDebug(className) ) return;
  	
  writeMachineInformation();
  writeTimeStamp();
  writeLogMessageType("debug");
  writeTraceInformation(className,methodName);

  _out << message << std::endl;
  _out.flush();
  }
}

void utils::CommandLineLogger::info(const std::string& className, const std::string& methodName, const std::string& message) {
  if ( WRITE_LOG_INFORMATION ) {
  writeMachineInformation();
  writeTimeStamp();
  writeLogMessageType("info");
  writeTraceInformation(className,methodName);

  _out << message << std::endl;
  _out.flush();
  }
}

void utils::CommandLineLogger::warning(const std::string& className, const std::string& methodName, const std::string& message) {
  if ( WRITE_LOG_INFORMATION ) {
  writeMachineInformation();
  writeTimeStamp();
  writeLogMessageType("warning");
  writeTraceInformation(className,methodName);

  _out << message << std::endl;
  _out.flush();
  }
}

void utils::CommandLineLogger::error(const std::string& className, const std::string& methodName, const std::string& message) {
  if ( WRITE_LOG_INFORMATION ) {
  writeMachineInformation();
  writeTimeStamp();
  writeLogMessageType("error");
  writeTraceInformation(className,methodName);

  _out << message << std::endl;
  _out.flush();
  }
}

void utils::CommandLineLogger::writeLogMessageType( const std::string& type ) {
  if (WRITE_LOG_MESSAGE_TYPE) {
    _out << type << "\t";
  }
}


void utils::CommandLineLogger::writeTraceInformation( const std::string& className, const std::string& methodName ) {
  if (WRITE_TRACE_INFORMATION) {
    _out << className << "::" << methodName << "\t";
  }
}


void utils::CommandLineLogger::writeTimeStamp() {
  if (WRITE_TIME_STAMP) {
    // calender time: create struct and get time from system
    time_t* timeStamp = new time_t();
    assert( timeStamp!=NULL );
    time(timeStamp);

    // Break down time into hour, seconds, ...
    // Note that time is only a substructure of timeStamp. Therefore the pointer
    // to time may not be deleted.
    tm*     time      = localtime(timeStamp);
    assert( time!=NULL );

    // write all information
    _out << time->tm_hour << ":" << time->tm_min << ":" << time->tm_sec << "\t";

    delete timeStamp;
  }
}

void utils::CommandLineLogger::writeMachineInformation() {
  if (WRITE_MACHINE_INFORMATION) {
    utsname* utsdata = new utsname();
    assert( utsdata!=NULL );
    uname(utsdata);

    _out << "[" //utsdata->sysname << ", "
         << utsdata->nodename;
  
    #ifdef Parallel
    _out << ",rank:" << parallel::Node::getInstance().getRank();
    #endif       
    
    _out << "]\t";
         
//         << utsdata->release << "]\t";

    delete utsdata;
  }
}
