#ifndef CGWASDATAIO_CLASS
#define CGWASDATAIO_CLASS

#include "globals.h"
#include "gwas/CGWASData.h"
#include "gwas/CMetaGWAS.h"

#include <string>

class CGWASDataIO {
	public:
		static void writeSummaryOutput(std::string const&, GWASData const&, GWASResults const&);
		static void writeFilteredPlinkFile(std::string const&, GWASData const&);
		static GWASResults readGWASResults(std::string const&);
		static void writeMetaResultsFile(std::string const&, CMetaResults const&);
};

#endif //CGWASDATAIO_CLASS
