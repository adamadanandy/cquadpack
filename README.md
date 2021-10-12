# CQUADPACK

This package includes the *QUADPACK* Fortran codes rewritten in C.
Unlike ports created by translators, this port reimplements the algorithms
to take advantage of C progralm structure and dynamic memory. Every effort
has been made to leave the low level code intact, except where some
opportunity to convert from unstructured blocks to more structured form
was taken.

Some changes to the calling conventions for various routines have also
been made. In particular, the low level integrators, such as **'dqk??.c'**,
are treated as functions returning a value. In the original code, these
integrators are subroutines passing the result back through a parameter.

QUADPACK is well documented in the publication, ``QUADPACK, A Subroutine
Package for Automatic Integration,'' by R. Piessens, et al., Springer-Verlag,
1980.

This port was initially done by ESSS/cquadpack, and C. Bond (http://www.crbond.com/)
This version of CQUADPACK is a perfectly replica of the QUADPACK with modifications
based on ESSS/cquadpack, which allows the workspace restored only in stack.

## DEPENDENCIES

dqc25o.c : LAPACK: dgtsv

## FILES

The files are distributed among several directories for convenience in
identification. For application development, a different directory
structure would be necessary. Feel free to move the files according to
your requirements.

### SOURCE FILES

The complete set of source and header files are in the directory: **src/**.

If you have access to the QUADPACK documentation, you will be able to
match files in this directory to the original FORTRAN files.

## CAVEATS

The author makes no warranty about the correctness of any translation. This
is just as well, because there is no guarantee that the original codes were
error free in the first place. Every effort has been made to assure that the
programs are functional and that they perform as expected.
