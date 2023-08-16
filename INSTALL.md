# basis installation guide #

Mark A. Caprio
Department of Physics, University of Notre Dame

+ 08/16/2023 (mac): Created.

----------------------------------------------------------------

This package is usually used as a submodule of a larger project.  However, if
you are doing development work on `basis`, you may want to set up standalone
compilation.

You will first need to clone the repositories for the ND makefile and for
various supporting modules.  Change to the intended parent directory for all the
repositories (including `basis`), e.g.:
  
  ~~~~~~~~~~~~~~~~
  % cd ~/code
  ~~~~~~~~~~~~~~~~

Then clone all the respositories (including `basis`, if you haven't already
cloned it):

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git clone https://github.com/nd-nuclear-theory/am.git
  % git clone https://github.com/nd-nuclear-theory/basis.git
  % git clone https://github.com/nd-nuclear-theory/fmt.git --branch submodule
  % git clone https://github.com/nd-nuclear-theory/mcutils.git
  % git clone https://github.com/nd-nuclear-theory/ndmakefile.git
  % git clone https://github.com/nd-nuclear-theory/ndconfig.git
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Compilation also requires the Boost, GSL, and Eigen libraries. See the
`INSTALL.md` file under the `ndconfig` repository for notes on installing or
accessing these libraries.

Then, provide symlinks to the make file and to the appropriate config file, e.g.:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % ln -s ../ndmakefile/makefile
  % ln -s ../ndconfig/config-gnu.mk config.mk
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See `ndconfig`'s `INSTALL.md` for further guidance on picking the config file.

To compile the various test codes:

  ~~~~~~~~~~~~~~~~
  % make programs-test
  ~~~~~~~~~~~~~~~~
