#ifndef CGWASDATAIO_CLASS
#define CGWASDATAIO_CLASS

#include "gin/globals.h"
#include "gin/gwas/CGWASData.h"

#include <string>

class CGWASDataIO {
	public:
		static void writeSummaryOutput(std::string const&, GWASData const&, GWASResults const&);
		static void writeFilteredPlinkFile(std::string const&, GWASData const&);
		static GWASResults readGWASResults(std::string const&);
};

#endif //CGWASDATAIO_CLASS
