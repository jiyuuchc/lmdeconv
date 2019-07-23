mex -setup c++
mex -R2018a -v CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' lmsamplermex.cpp
mex -R2018a -v CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' lmdeconvmex.cpp
mex -R2018a -v CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' lmdeconv3dmex.cpp
mex -R2018a -v CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' lmsampler3dmex.cpp
