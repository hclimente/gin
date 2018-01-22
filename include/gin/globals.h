#ifndef GLOBALS
#define GLOBALS

#include "types.h"

#define SKAT    0
#define CHI2    1
#define TREND   2

#define CONSISTENCY 0
#define BIC         1
#define AIC         2
#define AICc        3
#define mBIC        4

#define RED	    "\e[0;31m"
#define BLACK	"\e[0m"
#define GREEN	"\e[0;32m"
#define BLUE	"\e[0;34m"
#define YELLOW	"\e[0;33m"

#define LOG 	        ""
#define GIN_FERROR	    "FATAL ERROR"
#define GIN_ERROR	    "ERROR"
#define GIN_INFO	    "INFORMATION"
#define GIN_WARNING     "WARNING"
#define GIN_DEBUG	    "DEBUG"
#define GIN_STATUS	    "STATUS"
#define GIN_ATTENTION	"ATTENTION"

#define GIN_PI 3.14159265359f

//Writing log files
#ifndef AS_RGINLIB
#include <iostream>
#define logging(B,C) {\
	time_t rt; struct tm* ct;\
	time(&rt);\
	ct = localtime(&rt);\
	if((std::string)B==GIN_ERROR || (std::string)B==GIN_FERROR){\
		std::cerr << RED << "[" << ct->tm_mday << "." << ct->tm_mon+1 << "." << ct->tm_year + 1900\
		     << "," << ct->tm_hour << ":" << ct->tm_min << ":" << ct->tm_sec << "] " << B << " in "\
		     << __FILE__ << " at line " << __LINE__  << ": "\
		     << C << BLACK << "\n";\
	} else if((std::string)B==GIN_WARNING) {\
		std::cerr << YELLOW <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_INFO) {\
		std::cout << GREEN <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_STATUS) {\
		std::cout << BLUE <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_ATTENTION) {\
		std::cout << RED <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_DEBUG) {\
		std::cout << RED <<  C << BLACK << "\n";\
	} else {\
		std::cout << C << "\n";\
	}\
}
#else
#include "Rcpp.h"
#define logging(B,C) {\
	time_t rt; struct tm* ct;\
	time(&rt);\
	ct = localtime(&rt);\
	if((std::string)B==GIN_ERROR || (std::string)B==GIN_FERROR){\
		Rcpp::Rcerr << RED << "[" << ct->tm_mday << "." << ct->tm_mon+1 << "." << ct->tm_year + 1900\
		     << "," << ct->tm_hour << ":" << ct->tm_min << ":" << ct->tm_sec << "] " << B << " in "\
		     << __FILE__ << " at line " << __LINE__  << ": "\
		     << C << BLACK << "\n";\
	} else if((std::string)B==GIN_WARNING) {\
		Rcpp::Rcerr << YELLOW <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_INFO) {\
		Rcpp::Rcout << GREEN <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_STATUS) {\
		Rcpp::Rcout << BLUE <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_ATTENTION) {\
		Rcpp::Rcout << RED <<  C << BLACK << "\n";\
	} else if((std::string)B==GIN_DEBUG) {\
		Rcpp::Rcout << RED <<  C << BLACK << "\n";\
	} else {\
		Rcpp::Rcout << C << "\n";\
	}\
}
#endif //AS_RGINLIB

//Continued fraction struct
typedef struct cf_param {
	float64 a_odd;
	float64 b_odd;
	float64 a_even;
	float64 b_even;
} cf_param;

#endif //Globals