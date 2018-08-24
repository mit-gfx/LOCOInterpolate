# LOCOInterpolate
Algorithm for **LO**cal **CO**ntinuous Interpolation developed by Schulz et al 2017. 


# Setup
Clone repository:
```
git clone git@github.com:czwmit/schulzInterpolate.git
```

Download Eigen v3.3.5 from http://eigen.tuxfamily.org/ and place all Eigen files in this location:
```
LOCOInterpolate/libs/eigen/
```

Download CMake from https://cmake.org/

# Build on Unix/OSX

Download Google protobuf v3.1.0: https://github.com/google/protobuf/blob/master/src/README.md

or

Download Google protobuf v3.1.0 from: https://github.com/google/protobuf/releases/v3.1.0

Build protobuf following the instructions here: https://github.com/google/protobuf/blob/master/src/README.md

Place the libprotobuf.a file in 
```
LOCOInterpolate/libs/protobuf/lib
```

To build InstantCAD, execute the following in the schulzInterpolate directory:
```
$ mkdir build && cd build
$ cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release ..
$ make
$ cd ../../bin
& ./Test #To verify that the code works.
```

# Build on Windows using VisualStudio

Download and install Google protobuf v3.1.0 according to these directions https://github.com/google/protobuf/blob/master/cmake/README.md

**NOTE** - Build protobuf using a Win64 generator.

**NOTE** - Set the runtime library to Multi-threaded DLL /MD. In VisualStudio this can be done by Properties -> Configuration Properties -> C/C++ -> Code Generation -> Runtime Library -> Multi-threaded DLL /MD. If building in debug mode, then set the runtime library to Multi-threaded DLL debug /MDd.

Place the libprotobuf.lib file in 
```
LOCOInterpolate/libs/protobuf/lib
```

To build InstantCAD, execute the following in the schulzInterpolate directory:
```
C:\> mkdir build
C:\> cd build
C:\> cmake -G "Visual Studio [VERSION] Win64" -DCMAKE_BUILD_TYPE=Release ..
```

Open InstantCAD.sln and build the solution in x64 Release mode. Run Test.

# LCAdaptiveSampling

1. *OCAdaptiveGrid* represents the k-d tree in the parameter space. The tree stores samples, which are the result of evaluating various points in the parameter space. These samples are used to interpolate and approximate any given point in the space. Our algorithm adaptively refines this tree until the interpolation approximates the evaluation function to a desired degree of accuracy.
2. *OCAdaptiveGridCell* represents an element of the k-d tree
3. *OCBasisFunction* represents the basis function that compose the support of a given sample. This class is subclassed by LCLinearBSpline and LCCubicBSpline. These subclasses implement linear and cubic supports respectively. 
4. *OCSampledFunction* inherits from OCFunction and can evaluate any point in its domain. Note that this is accomplished using interpolation and not direct evaluation of the function.
5. *OCSample* represents the result of evaluating a particular coordinate in the parameter space.
6. *OCProtoConverter* handles storage of samples and functions via Google's protobuf 

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

[^fn]: http://people.csail.mit.edu/aschulz/instantCAD/index.html
