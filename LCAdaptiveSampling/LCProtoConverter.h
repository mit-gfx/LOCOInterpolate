#ifndef OC_PROTO_CONVERTER_H
#define OC_PROTO_CONVERTER_H
#include "LCShapeInfoProtoConverter.h"
#include "PrecomputedParametricShape.pb.h"
class LCSampledFunction;
using namespace OptCAD;

class LCError;


class PrecomputedParametricShapeConverter : public ProtoConverter<proto::PrecomputedParametricShape>
{
public:

	static LCError outputToProto(std::string filename, bool useAscii, LCSampledFunction* shape);
	static LCError loadFromProto(std::string filename, bool useAscii, LCSampledFunction** shape);
private:	
	static LCError convertToProto(LCSampledFunction* paramShape, std::string filepath, proto::PrecomputedParametricShape* result);
	static LCError convertToCplusplus(const proto::PrecomputedParametricShape& protoParamShape,
		LCSampledFunction** result);
};


#endif
