%module CEasyGWAS

/*
*Include all headers
*/
%{
#define SWIG_FILE_WITH_INIT
#define SWIG
#include "globals.h"
#include "types.h"
#include "regression/CRegression.h"
#include "gwas/CSingleTraitGWAS.h"
#include "gwas/CGWASData.h"
#include "gwas/CMetaGWAS.h"
#include "meta/CMetaAnalysis.h"
#include "meta/CEffectSize.h"
#include "kernel/CKernels.h"
#include "stats/CChi2.h"
#include "stats/CGaussian.h"
#include "stats/CBeta.h"
#include "stats/CFisherF.h"
#include "stats/CStudentT.h"
#include "stats/CGamma.h"
#include "stats/CHSIC.h"
/*
*/
%}

%include "types.i"

/*
*Include std lib specific swig interfaces
*/
%include "exception.i"
%include "std_string.i"
%include "std_map.i"
%include "std_vector.i"

namespace std {
        %template(StringVector) vector<string>;
        %template(CharVector) vector<char>;
        %template(CharCharVector) vector< vector<char> >;
        %template(Uint64Vector) vector<uint64>;
        %template(UintVector) vector<uint>;
        %template(Float64Vector) vector<float64>;
}

/*
*Add external numpy interfaces to access python specific numpy types
*/
%include "External/numpy.i"


%init %{
        import_array();
%}

/*
*Include Eigen lib to convert between numpy and eigen types
*/
%include "eigen.i"


/*
*Include internal interfaces
*/
%include "regression/CRegression.i"
%include "gwas/CSingleTraitGWAS.i"
%include "meta/CMetaAnalysis.i"
%include "meta/CEffectSize.i"
%include "gwas/CGWASData.i"
%include "gwas/CMetaGWAS.i"
%include "kernel/CKernels.i"
%include "stats/CChi2.i"
%include "stats/CGaussian.i"
%include "stats/CBeta.i"
%include "stats/CFisherF.i"
%include "stats/CStudentT.i"
%include "stats/CGamma.i"
%include "stats/CHSIC.i"
/*
*Ignore some global stuff
*/
