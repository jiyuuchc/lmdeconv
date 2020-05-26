## lmdeconv
The software implements a few algorithms for single molecule localization microscopy. It includes a sampler, a function to perform MAP estiamtion of molecular distribution, and a particle fusion algorithm.

**Ref:** 10.1073/pnas.1912634116

### Requirements
Matlab with image processing and optimization toolbox. I've only tested on Linux (Debian) and Windows 10 (64bit).

### Usage: Paricle Fusion workflow

##### 1. Build
Make sure you have gcc C-comiler installed. On windows, it is easiest to use the "AddOn..." menu option, and search for Cygwin. 
In Matlab, run
```
build
```
##### 2. Load data. 
Read your particle data into a cell array (e.g., particles). Each cell stores one particle using a Nx3 array. Each row of array is (x, y, sigma). Sigma is an estimate of localization error. Use consitent unit for all values. 

##### 3. Pre-register against one particle (optional)
```
target = particle{1} % choose any particle, presumbly one with good data quality 
pixelsize = 1.0 % choose a pixelsize for rendering. Same unit as your particle data
nsamples = 5000 % how many samples to be drawn to build a alignement template
particles = particleIter(particles, pixelsize, nsamples, 'Target', target)
```

##### 4. Register iteration
Performing global registration.
```
particles = particleIter(particles, pixelsize, nsamples)
```
This command should be run multiple times, which improve alignment iteratively.

##### 5. Check fusion results
```
dataobj = lmdatainit(cat(1,particles{:}), pixelsize)
imagesc(histimg(dataobj))
axis off equal
```

##### 6. Perform MAP rendering of the fusion result
```
imagesc(lmdeconv(dataobj, 200))
axis off equal
```
