<h1 align="center">libff</h1>
<h4 align="center">C++ library for Finite Fields and Elliptic Curves</h4>

___libff___ is a C++ library for finite fields and elliptic curves. The library is developed by [SCIPR Lab] and contributors (see [AUTHORS] file) and is released under the MIT License (see [LICENSE] file).

## Table of contents

- [Directory structure](#directory-structure)
- [Elliptic curve choices](#elliptic-curve-choices)
- [Build guide](#build-guide)

## Directory structure

The directory structure is as follows:

* [__libff__](libff): C++ source code, containing the following modules:
  * [__algebra__](libff/algebra): fields and elliptic curve groups
  * [__common__](libff/common): miscellaneous utilities
* [__depends__](depends): dependency libraries

## Elliptic curve choices

The libsnark library currently provides three options:

* `edwards`:
   an instantiation based on an Edwards curve, providing 80 bits of security.

* `bn128`:
   an instantiation based on a Barreto-Naehrig curve, providing 128
   bits of security. The underlying curve implementation is
   \[ate-pairing], which has incorporated our patch that changes the
   BN curve to one suitable for SNARK applications.

    *   This implementation uses dynamically-generated machine code for the curve
        arithmetic. Some modern systems disallow execution of code on the heap, and
        will thus block this implementation.

        For example, on Fedora 20 at its default settings, you will get the error
        `zmInit ERR:can't protect` when running this code. To solve this,
        run `sudo setsebool -P allow_execheap 1` to allow execution,
        or use `make CURVE=ALT_BN128` instead.

* `alt_bn128`:
   an alternative to `bn128`, somewhat slower but avoids dynamic code generation.

Note that `bn128` requires an x86-64 CPU while the other curve choices
should be architecture-independent.

## Build guide

The library has the following dependencies:

* [Boost](http://www.boost.org/)
* [CMake](http://cmake.org/)
* [GMP](http://gmplib.org/)
* [libprocps](http://packages.ubuntu.com/trusty/libprocps-dev)

The library has been tested on Linux, but it is compatible with Windows and Mac OS X.

### Installation

On Ubuntu 14.04 LTS:

```
sudo apt-get install build-essential git libboost-all-dev cmake libgmp3-dev libssl-dev libprocps3-dev pkg-config
```

Fetch dependencies from their GitHub repos:

```
git submodule init && git submodule update
```

### Compilation

To compile, starting at the project root directory, create the build directory and Makefile:

```
mkdir build && cd build && cmake ..
```
Optionally, you can specify the install location by providing the desired install path prefix:
```
cmake .. -DCMAKE_INSTALL_PREFIX=/install/path
```

Then, to compile and install the library, run this within the build directory:
```
make
make install
```

This will install `libff.a` into `/install/path/lib`; so your application should be linked using `-L/install/path/lib -lff`. It also installs the requisite headers into `/install/path/include`; so your application should be compiled using `-I/install/path/include`.

### Bazel

```
$ bazel query 'attr(visibility, "//visibility:public", libff/...:*)'
$ bazel build libff
$ bazel build libff/algebra/...
...etc...
$ bazel test test
```

To use in an app, add the following to your WORKSPACE file.

```
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
http_archive(
    name = "rules_foreign_cc",
    strip_prefix="rules_foreign_cc-master",
    url = "https://github.com/bazelbuild/rules_foreign_cc/archive/master.zip",
)
load("@rules_foreign_cc//:workspace_definitions.bzl",
     "rules_foreign_cc_dependencies")
rules_foreign_cc_dependencies()

## Bazelized external repos:
local_repository( name = "libff" , path = "depends/libff")
# http_archive(
#     name = "libff",
#     urls = ["https://github.com/obazl/libff/archive/bzl-1.0.tar.gz"],
#     strip_prefix = "libff-bzl-1.0",
#     sha256 = ...obtain from github...
# )
local_repository( name = "ate_pairing" , path = "depends/ate-pairing")
# http_archive(
#     name = "ate_pairing",
#     urls = ["https://github.com/obazl/ate-pairing/archive/bzl-1.0.tar.gz"],
#     strip_prefix = "ate-pairing-bzl-1.0",
#     sha256 = ...obtain from github...
# )
## build target: @xbyak//xbyak
local_repository( name = "xbyak" , path = "depends/xbyak")
# http_archive(
#     name = "xbyak",
#     urls = ["https://github.com/obazl/xbyak/archive/bzl-1.0.tar.gz"],
#     strip_prefix = "xbyak-bzl-1.0",
#     sha256 = ...obtain from github...
# )
## Non-bazelized external repos:
all_content = """filegroup(name = "all", srcs = glob(["**"]), visibility = ["//visibility:public"])"""
http_archive(
    name="openmp",
    url="https://github.com/llvm/llvm-project/releases/download/llvmorg-10.0.0/openmp-10.0.0.src.tar.xz",
    sha256="3b9ff29a45d0509a1e9667a0feb43538ef402ea8cfc7df3758a01f20df08adfa",
    strip_prefix="openmp-10.0.0.src",
    build_file_content = all_content
)
http_archive(
    name="openssl",
    url="https://www.openssl.org/source/openssl-1.1.1g.tar.gz",
    sha256="ddb04774f1e32f0c49751e21b67216ac87852ceb056b75209af2443400636d46",
    strip_prefix="openssl-1.1.1g",
    build_file_content = all_content
)
http_archive(
    name="libgmp",
    url="https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz",
    sha256="258e6cd51b3fbdfc185c716d55f82c08aff57df0c6fbd143cf6ed561267a1526",
    strip_prefix = "gmp-6.2.0",
    build_file_content = all_content
)
```

## Testing

To execute the tests for this library, run:
```
make check
```

## Profile

To compile the multi-exponentiation profiler in this library, run:
```
make profile
```
The resulting profiler is named `multiexp_profile` and can be found in the `libff` folder under the build directory.

[SCIPR Lab]: http://www.scipr-lab.org/ (Succinct Computational Integrity and Privacy Research Lab)

[LICENSE]: LICENSE (LICENSE file in top directory of libff distribution)

[AUTHORS]: AUTHORS (AUTHORS file in top directory of libff distribution)
