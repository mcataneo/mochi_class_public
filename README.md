mochi_class: Modelling Optimisation to Compute Horndeski in CLASS
==============================================

<p align="left">
  <img src="logos/mochi_class_logo_github_black_back.png" alt="mochi_class logo" width="300" />
</p>

Authors: Matteo Cataneo, Emilio Bellini

This code is an extension of the hi_class code (http://miguelzuma.github.io/hi_class_public/), itself a patch to the Einstein-Boltzmann solver CLASS (http://class-code.net). The new features include:

  * $\alpha$-functions replaced by stable basis {$M_{\m Pl}^2$, $D_{\rm kin}$, $c_{\rm s}^2$}
  *
  *
  *


Compiling mochi_class and getting started
-----------------------------------

Clone the code from [https://github.com/mcataneo/mochi_class_public](https://github.com/mcataneo/mochi_class_public.git). 
Go to the mochi_class directory (cd mochi_class_public/) and compile (make clean;
make class). You can usually speed up compilation with the option -j:
make -j class. If the first compilation attempt fails, you may need to
open the Makefile and adapt the name of the compiler (default: gcc-13),
of the optimization flag (default: -O3) and of the OpenMP
flag (default: -fopenmp; this flag is facultative, you are free to
compile without OpenMP if you don't want parallel execution). 
Many more details (although a bit out-of-date) on the CLASS compilation are given on the
wiki page

https://github.com/lesgourg/class_public/wiki/Installation

(in particular, for compiling on Mac >= 10.9 despite of the clang
incompatibility with OpenMP).

To test the code has been properly installed, type, for instance:

    ./class /inifiles/galileon.ini

If everything works as expected, this will generate CMB and matter power spectra in the /output directory.

The `<cosmology>.ini` files in the inifiles/ directory can be used as reference input files, containing and
explaining the use of all possible input parameters. When creating
your own input files, make sure they have a *.ini
extension.

If you want to play with the precision/speed of the code, you can use
one of the provided precision files (e.g. cl_permille.pre) or modify
one of them, and run with two input files, for instance:

    ./class test.ini cl_permille.pre

The files *.pre are suppposed to specify the precision parameters for
which you don't want to keep default values. If you find it more
convenient, you can pass some (or all of) these precision parameter values in your *.ini
file instead of an additional *.pre file.

Python
------

To use mochi_class from python, or ipython notebooks, or from the Monte
Python parameter extraction code, you need to compile not only the
code, but also its python wrapper. This can be done by typing just
'make' instead of 'make class' (or for speeding up: 'make -j'). To ensure 
a smooth experience, before installing class we recommend you first install 
its dependencies:

  * numpy <= 1.26.4
  * cython >= 3.x.y
  * setuptools <= 65.5.1

We provide an example notebook (stable_params.ipynb) using the new mochi_class features. It can be found in the /notebooks directory.

More details on the wrapper and its compilation are found on the wiki page (again slightly out-of-date)

https://github.com/lesgourg/class_public/wiki


Using the code
--------------

You can use CLASS freely, provided that in your publications, you cite
at least the paper `CLASS II: Approximation schemes <http://arxiv.org/abs/1104.2933>`. Feel free to cite more CLASS papers!

Support
-------

To get support, please open a new issue on the

https://github.com/lesgourg/class_public

webpage!
