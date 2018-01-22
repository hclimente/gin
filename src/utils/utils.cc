//
// Created by hclimente on 06/10/2017.
//

#include "gin/utils/utils.h"

int libgin_present() {
	return 1;
}

#ifndef AS_RGINLIB
#include <stdlib.h>
void abort_gin(int status) {
	exit(status);
}
#else
#include "Rcpp.h"
void abort_gin(int status) {
	Rcpp::stop(status);
}
#endif //AS_RGINLIB
