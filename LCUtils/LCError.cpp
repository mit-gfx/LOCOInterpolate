#include "LCError.h"

LCError::LCError()
{
	isOK_ = true;
	internalDescription_ = "no error";
}
LCError::LCError(std::string internalDescription)
{
	isOK_ = false;
	internalDescription_ = internalDescription;
}
std::string LCError::internalDescription() const
{
	return internalDescription_;
}
bool LCError::isOK() const
{
	return isOK_;
}
