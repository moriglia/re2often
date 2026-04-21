(use-modules (guix-science packages fortran)
             (gnu packages commencement))

(packages->manifest
 (list fortran-fpm
       gfortran-toolchain))
