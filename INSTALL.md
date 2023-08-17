# basis installation guide #

Mark A. Caprio  
Department of Physics, University of Notre Dame  
Patrick J. Fasano  
Physics Division, Argonne National Laboratory

+ 08/16/2023 (mac): Created.
+ 08/17/2023 (pjf): Switch to CMake.

----------------------------------------------------------------

This package is usually used as a submodule of a larger project.  However, if
you are doing development work on `basis`, you may want to set up standalone
compilation.

You will first need to build (and optionally install) the various supporting
modules.  Change to the intended parent directory for all the repositories
(including `basis`), e.g.:

  ~~~~~~~~~~~~~~~~
  % cd ~/code
  ~~~~~~~~~~~~~~~~

Then clone all the repositories (including `basis`, if you haven't already
cloned it):

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git clone https://github.com/nd-nuclear-theory/am.git
  % git clone https://github.com/nd-nuclear-theory/mcutils.git
  % git clone https://github.com/fmtlib/fmt.git
  % git clone https://github.com/nd-nuclear-theory/basis.git
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each of `am`, `mcutils`, and `fmt`, you should use CMake to build (and
optionally install) the library.

Compilation also requires the Boost, GSL, and Eigen libraries. See the
`INSTALL.md` file under the `ndconfig` repository for notes on installing or
configuring access to these libraries.

Next, configure the project using CMake using:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cmake -B build/ .
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have not installed the dependencies, but, e.g., only built them, you may
need to tell CMake where to find the built-but-not-installed libraries:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cmake -B build . \
      -Dam_DIR=~/code/am/build \
      -Dmcutils_DIR=~/code/mcutils/build \
      -Dfmt_DIR=~/code/fmt/build
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compile the library itself:

  ~~~~~~~~~~~~~~~~
  % cmake --build build/
  ~~~~~~~~~~~~~~~~

To compile the various test codes:

  ~~~~~~~~~~~~~~~~
  % cmake --build build/ -- tests
  ~~~~~~~~~~~~~~~~

To install the library (here with prefix `~/install`):
  ~~~~~~~~~~~~~~~~
  % cmake --install build/ --prefix ~/install
  ~~~~~~~~~~~~~~~~
