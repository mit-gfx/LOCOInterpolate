#ifndef OC_SHAPEINFO_PROTO_CONVERTER_H
#define OC_SHAPEINFO_PROTO_CONVERTER_H
#include "LCError.h"
#include <google/protobuf/text_format.h>

#include "PrecomputedParametricShape.pb.h"
using namespace OptCAD;

class LCError;

template<class ProtoClass>
class ProtoConverter
{
protected:
	static LCError saveToFile(const ProtoClass& proto, std::string filename, bool ascii = false)
	{
		LCError err;
		if (ascii) {
			std::string output;
			bool success = google::protobuf::TextFormat::PrintToString(proto, &output);
			if (!success)
			{
				return LCError("could not print to proto");
			}
			std::ofstream outfile(filename, std::ios::out);
			if (!outfile)
			{
				return LCError("could not write to file " + filename);
			}
			outfile << output;
			outfile.close();
		}
		else {
			std::ofstream outfile(filename, std::ios::out | std::ios::binary);
			if (!outfile)
			{
				return LCError("could not write to file " + filename);
			}
			proto.SerializePartialToOstream(&outfile);
			outfile.close();
		}
		return err;
	}

	
	static LCError loadFromFile(std::string filename, bool ascii, ProtoClass* result)
	{		
		LCError err;
		if (ascii) {
			std::ostringstream input;
			std::ifstream infile(filename, std::ios::in);
			if (!infile)
			{
				return LCError("could not open file " + filename);
			}
			input << infile.rdbuf();
			std::string inputStr = input.str();
			bool success = google::protobuf::TextFormat::ParseFromString(inputStr, result);
			if (!success)
			{
				return LCError("could not parse the file");
			}
		}
		else {
			std::ifstream infile(filename, std::ios::in | std::ios::binary);
			if (!infile)
			{
				return LCError("could not open file " + filename);
			}
			bool success = result->ParseFromIstream(&infile);
			if (!success)
			{
				return LCError("could not parse the file");
			}
		}
		return err;
	}


};


/*removed by czw
class OCTetCADShapeInfoConverter : public ProtoConverter<proto::TetMesh>
{	
public:

	static LCError outputToProto(std::string filename, OCTetCADShapeInfo* shape);
	static LCError loadFromProto(std::string filename, OCTetCADShapeInfo* shape);
private:
	static LCError convertToProto(OCTetCADShapeInfo* paramShape, proto::TetMesh* result);
	static LCError convertToCplusplus(const proto::TetMesh & protoParamShape,
		OCTetCADShapeInfo* result);
};
*/

#endif
