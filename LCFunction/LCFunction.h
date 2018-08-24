#ifndef _OC_PARAMETRIC_SHAPE_
#define _OC_PARAMETRIC_SHAPE_


#include <vector>
#include <Eigen/Dense>

class LCError;
class LCFunctionValue;
class LCFunctionDeriv;

/*  This is an abstract class that defines a parametric shape
	A parametric shape is abstractly defined by a feasible set and a function that maps point
	in this feasible set to shapes. In our code, the feasible set is always a hypercube defined by a set
	of ranges. The functions below return information about the feasible set (validate parameters, help sample
	the feasible set, etc) and also return the shape (and the derivative of the shape) for a given 
	parameter's configuration.
 */
class LCFunction {
public:
	LCFunction();

	/* returns the number of parameters */
	virtual int getNParams() const = 0;
	/* gets the min range for param i - does not check existence */
	virtual double getMinRange(int iParam) const = 0;
	/* gets the max range for param i - does not check existence */
	virtual double getMaxRange(int iParam) const = 0; 
	/*returns the shape for a give set of parameters*/
	virtual LCError evalShapeInfo(const std::vector<double> &params, LCFunctionValue **result) = 0;
	/*returns the shape derivative for a give set of parameters*/
	virtual LCError evalDeriv(const std::vector<double> &params, int direction, LCFunctionDeriv **result) = 0;


	/* returns the parameters at the center of the paremeter space */
	LCError getMidPoint(std::vector<double> *result);
	/* checks if parameters are valid */
	LCError validateParameters(const std::vector<double> &parameters) const;
	/* maps a value between 0 and 1 to the ranges of param i */
	LCError getRelativeParamVal(int i, double alpha, double *retult);
	/* maps a point in the hypercube with ranges [-1 1] to the parameter space (given by the ranges) */
	std::vector<double> mapToStandardHypercube(const std::vector<double> &params);
	/* maps a point in paramter space to the hypercube with ranges [-1 1] */
	std::vector<double> mapFromStandardHypercube(const std::vector<double> &unitParams);
	/* returns the vector of ranges */
	void getRanges(std::vector<double> * ranges);

};

#endif
