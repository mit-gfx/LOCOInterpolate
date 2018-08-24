# LCFunction

This library contains the abstract class that defines a a function _f_: _X_-> _Y_, where:

* _X_ is the parameter space and is contained in _R_^_N_ (_N_= number of parameters). More precisely, _X_ is an _N_-dimensional rectangular prism with boundaries given by `LCFunction.getRanges()`.

* _Y_ is the image space, and can represent a variety of abstract values (e.g. it can be a real number or a 3D mesh). This domain is not implemented in the base class, but will be defined by the subclasses which inherit from OCFunction.

The four base classes are:

* *LCFunction* which represents a function _f_ on the parameter space.
* *LCFunctionValue* represents the result of evaluating _f_ at a particular coordinate
* *LCFunctionDeriv* represents the result of evaluating the derivative of _f_ at a particular coordinate
* *LCFunctionValueMap* specifies a homeomorphic mapping between two values in the image of _f_

We also provide a set of four classes that inherit from these base classes. Together, these classes define a function _f_ : _X_ -> _R_. Users who wish to define their own functions and image space should modify these files or follow their example:

* LCRealFunction
* LCRealFunctionValue
* LCRealFunctionDerivValue
* LCRealFunctionValueMap

Finally, the classes PrecomputedParametricShape.pb and LCShapeInfoProtoConverter handle data storage via Google's protobuf.
