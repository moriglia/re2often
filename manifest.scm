;; In order to use this manifest, make sure the guix-hpc channel
;; is included in your channels file.
;; If not, you can include the following in your `$HOME/.config/guix/channels.scm`
;; (append
;;  (list
;;   ;; --- Other channels you may already have
;;   ;; (channel ...)
;;   ;; (channel ...)
;;   ;;
;;   ;; --- Guix-HPC channel
;;   (channel
;;    (name 'guix-hpc)
;;    (url "https://gitlab.inria.fr/guix-hpc/guix-hpc.git")
;;    (branch "master")))
;;  %default-channels)


(use-modules (guix packages)
	     (guix git-download)
	     (guix build-system cmake)
	     (guix build-system copy)
	     (guix licenses)
	     (gnu packages)
	     (gnu packages gcc)
	     (gnu packages cmake)
	     (gnu packages mpi)
	     (guix gexp)
	     (guix download)
	     (guix profiles)
	     (guix-science packages fortran)
	     (guix-hpc packages toolchains))


(define-public libprimaf
  (package
   (name "libprimaf")
   (version "0.7.2")
   (source
    (origin
     (method git-fetch)
     (uri (git-reference
           (url "https://github.com/libprima/prima.git")
           (commit (string-append "v" version))))
     (sha256
      (base32
       "1mcw5g2sj84j52300c2v8dlrzi6gmdyzddhi7g7k7kxg2ihwp5r5"))))
   (build-system cmake-build-system)
   (arguments
    `(#:tests? #f
      #:configure-flags
      '("-DCMAKE_INSTALL_LIBDIR=lib"
	"-DCMAKE_INSTALL_INCLUDEDIR=include")))
   (inputs
    `(("gfortran-toolchain" ,gfortran)))
   (synopsis "PRIMA: Reference Implementation for Powell's Methods with Modernization and Amelioration")
   (description
    "PRIMA is a package for solving general nonlinear optimization problems without using derivatives. It provides the reference implementation for Powell's renowned derivative-free optimization methods, i.e., COBYLA, UOBYQA, NEWUOA, BOBYQA, and LINCOA. The \"P\" in the name stands for Powell, and \"RIMA\" is an acronym for \"Reference Implementation with Modernization and Amelioration\".")
   (home-page "https://www.libprima.net")
   (license bsd-3)))

(define-public open-coarrays
  (package
   (name "open-coarrays")
   (version "2.10.2")
   (source
    (origin
     (method git-fetch)
     (uri (git-reference
           (url "https://github.com/sourceryinstitute/OpenCoarrays.git")
           (commit version)))
     (sha256
      (base32
       "1qcbqmj1kpd8lhjva6sb7250rhjw4s6ly553gssic2zc6w1zzmg7"))))
   (build-system cmake-build-system)
   (arguments
    `(#:tests? #f))
   (inputs
    `(("gfortran-toolchain" ,gfortran-14)
      ("mpi" ,openmpi)))
   (synopsis "Open Coarrays library")
   (description
    "OpenCoarrays is an open-source software project that produces an application
binary interface (ABI) used by the GNU Compiler Collection (GCC) Fortran
front-end to build executable programs that leverage the parallel programming
features of Fortran 2018.")
   (home-page "https://www.opencoarrays.org")
   (license bsd-3)))



(packages->manifest
 (list gfortran-toolchain-14
       open-coarrays
       libprimaf
       fortran-fpm))
