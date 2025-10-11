sr-indexes: A Fast and Small Subsampled r-index and r-CSA
=====

Getting Started
-----

### Downloading
To download our project, execute:
```shell
$ git clone https://github.com/duscob/sr-index.git
```

### Compiling
sr-index is a C++ template library using standard C++14.
As a build system, the project uses [CMake](https://cmake.org/).

>##### Dependencies
>
>The sr-index library requires the following external libraries:
>  * [Succinct Data Structure Library (SDSL)](https://github.com/simongog/sdsl-lite "SDSL's GitHub repository")
>
>`SDSL` must be installed on your system.


##### Building

To build our solution, you can create a build folder on the source directory and move to it.
```shell
$ cd dret
$ mkdir build
$ cd build
```

Then build and compile the project using the commands `cmake` and `make`, respectively.
```shell
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
```
