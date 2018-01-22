//
// Created by hclimente on 06/10/2017.
//

#include "gin/utils/utils.h"
#include "gin/utils/StringHelper.h"

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
	Rcpp::stop(StringHelper::to_string<int>(status));
}
#endif //AS_RGINLIB
