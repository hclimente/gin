#include "gin/io/CLogging.h"
#include "gin/utils/utils.h"

CLogging::CLogging() {
	file_logging = false;
}

CLogging::CLogging(std::string const& filename) {
	__filename = filename;
	__ofs.open(filename.c_str());
	if(!__ofs.is_open()) {
		logging(GIN_ERROR,"Could not open log file. Only screen logging supported!");
		abort_gin(-1);
	} else {
		file_logging = true;
	}
}

CLogging::~CLogging() {
	if(file_logging)
		__ofs.close();
}

void CLogging::log(std::string const& mode, std::string const& msg) {
	if(file_logging) {
		logging(mode,msg);
		if(mode==GIN_ERROR) {
			time_t rt;
			struct tm* ct;
			time(&rt);
			ct = localtime(&rt);
			__ofs << "[" << ct->tm_mday << "." << ct->tm_mon+1 << "." << ct->tm_year + 1900\
		     		<< "," << ct->tm_hour << ":" << ct->tm_min << ":" << ct->tm_sec << "] " << mode << " in "\
		     		<< __FILE__ << " at line " << __LINE__  << ": "\
		     		<< msg << "\n";\
		} else {
			__ofs << msg << "\n";
		}	
	} else {
		logging(mode,msg);
	}
}
