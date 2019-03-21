mex -setup c++
mex -R2018a -v CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' lmsamplermex.cpp
mex -R2018a -v CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' lmdeconvmex.cpp
