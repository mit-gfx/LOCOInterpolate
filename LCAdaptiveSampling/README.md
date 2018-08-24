# LCAdaptiveSampling

1. *OCAdaptiveGrid* represents the k-d tree in the parameter space. The tree stores samples, which are the result of evaluating various points in the parameter space. These samples are used to interpolate and approximate any given point in the space. Our algorithm adaptively refines this tree until the interpolation approximates the evaluation function to a desired degree of accuracy.
2. *OCAdaptiveGridCell* represents an element of the k-d tree
3. *OCBasisFunction* represents the basis function that compose the support of a given sample. This class is subclassed by LCLinearBSpline and LCCubicBSpline. These subclasses implement linear and cubic supports respectively. 
4. *OCSampledFunction* inherits from OCFunction and can evaluate any point in its domain. Note that this is accomplished using interpolation and not direct evaluation of the function.
5. *OCSample* represents the result of evaluating a particular coordinate in the parameter space.
6. *OCProtoConverter* handles storage of samples and functions via Google's protobuf 
