#ifndef _OC_ERROR_
#define _OC_ERROR_
#define LCErrorReturn(err) { if (!err.isOK())	{ return err; }}

#include <string>

class LCError
{
public:

	LCError();
	LCError(std::string internalDescription);
	std::string internalDescription() const;
	bool isOK() const;


private:
	bool isOK_;
	std::string internalDescription_;
};

#endif