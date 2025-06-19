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
	     (guix-hpc packages toolchains))


(define-public fortran-fpm
  (package
   (name "fortran-fpm")
   (version "0.12.0")
   (source
    (origin
     (method git-fetch)
     (uri (git-reference
           (url "https://github.com/fortran-lang/fpm.git")
           (commit version)))
     (sha256
      (base32
       "1lvmf8w3wfqyvv1r4jsslq58j2b35gfffyz605vfsapkcr8jzwc0"))))
   (build-system copy-build-system)
   (arguments
    (list
     #:phases
     #~(modify-phases
	%standard-phases
	(add-after 'patch-generated-file-shebangs 'build
		   (lambda* (#:key #:allow-other-keys)
			    (let* ((source (assoc-ref %build-inputs "source"))
				   (build-dir (string-append source "/../build"))
				   (bootstrap-dir
				    (string-append build-dir "/bootstrap"))
				   (fpm-bootstrap-file
				    #$(this-package-input "bootstrap"))
				   (out (assoc-ref %outputs "out")))
			      (mkdir-p bootstrap-dir)
			      (invoke "gfortran"
				      "-J" bootstrap-dir
				      "-o" (string-append bootstrap-dir "/fpm")
				      fpm-bootstrap-file)
			      (invoke (string-append bootstrap-dir "/fpm")
				      "install"
				      "--prefix" out))))
	(add-before
	 'build 'patch-fpm-toml
	 (lambda* (#:key inputs #:allow-other-keys)
		  (use-modules (ice-9 textual-ports) (srfi srfi-1))
		  (let* ((deps-dir "local-deps")
			 (deps `(("toml-f" . "toml-f")
				 ("M_CLI2" . "M_CLI2")
				 ("fortran-regex" . "fortran-regex")
				 ("jonquil" . "jonquil")
				 ("fortran-shlex" . "fortran-shlex"))))
		    (mkdir-p deps-dir)
		    ;; Symlink each dependency into local-deps/
		    (for-each
		     (lambda (dep)
		       (let* ((name (car dep))
			      (dir (assoc-ref inputs name)))
			 (symlink dir (string-append deps-dir "/" (cdr dep)))))
		     deps)
		    (substitute*
		     "fpm.toml"
		     (("toml-f\\.git.*")
		      (string-append "toml-f.path = \"" deps-dir "/toml-f\"\n"))
		     (("M_CLI2\\.git.*")
		      (string-append "M_CLI2.path = \"" deps-dir "/M_CLI2\"\n"))
		     (("fortran-regex\\.git.*")
		      (string-append "fortran-regex.path = \"" deps-dir "/fortran-regex\"\n"))
		     (("jonquil\\.git.*")
		      (string-append "jonquil.path = \"" deps-dir "/jonquil\"\n"))
		     (("fortran-shlex\\.git.*")
		      (string-append "fortran-shlex.path = \"" deps-dir "/fortran-shlex\"\n"))
		     (("^.*\\.rev.*") "")
		     (("^.*\\.tag.*") "")
		     )))
	 )
	(delete 'install)
	)
     )
    )
   (inputs
    `(("gfortran-toolchain" ,gfortran)
      ("bootstrap"
       ,(origin
	 (method url-fetch)
	 (uri "https://github.com/fortran-lang/fpm/releases/download/v0.10.1/fpm-0.10.1.F90")
	 (sha256
	  (base32
	   "1jp1hh5yjz3ghk5f1ixz06x6jl8mhh3hm89vla4yi9y2c1dx0lvm"))
	 (file-name "fpm.F90")))
      ;; ("git" ,git)
      ("toml-f"
       ,(origin
	 (method git-fetch)
	 (uri (git-reference
	       (url "https://github.com/toml-f/toml-f")
	       (commit "d7b892b1d074b7cfc5d75c3e0eb36ebc1f7958c1")))
	 (sha256
	  (base32
	   "07r78dk0q7xxrh3fjfrsx5jmf6ap21b4d5qld6ga3awki0cm75z8"))))
      ("M_CLI2"
       ,(origin
	 (method git-fetch)
	 (uri (git-reference
	       (url "https://github.com/urbanjost/M_CLI2.git")
	       (commit "7264878cdb1baff7323cc48596d829ccfe7751b8")))
	 (sha256
	  (base32
	   "06y7zndb0qcnyq7rxwhq0vv0j4dzkdw9hqphkp13k1zqf3vz8z28"))))
      ("fortran-regex"
       ,(origin
	 (method git-fetch)
	 (uri (git-reference
	       (url "https://github.com/perazz/fortran-regex")
	       (commit "1.1.2")))
	 (sha256
	  (base32
	   "183vxa1082kkg48rl75nxkm8j67vxpak3347dfzfbbxi0wyfklba"))))
      ("jonquil"
       ,(origin
	 (method git-fetch)
	 (uri (git-reference
	       (url "https://github.com/toml-f/jonquil")
	       (commit "4fbd4cf34d577c0fd25e32667ee9e41bf231ece8")))
	 (sha256
	  (base32
	   "1zk3rpl5npk4qaslcbp9nay6p9dsl3sqv2m0zc6337mxqa4dmjm0"))))
      ("fortran-shlex"
       ,(origin
	 (method git-fetch)
	 (uri (git-reference
	       (url "https://github.com/perazz/fortran-shlex")
	       (commit "2.0.0")))
	 (sha256
	  (base32
	   "0sygyjqwxyh974bmp8n5bxjs9nsdfavy5k3vhmx2v9bbl31jqk8a"))))
      ))
   (synopsis "Fortran Package Manager")
   (description
    "Fortran Package Manager (fpm) is a package manager and build system for Fortran.
Its key goal is to improve the user experience of Fortran programmers.
It does so by making it easier to build your Fortran program or library,
run the executables, tests, and examples, and distribute it as a dependency to other Fortran projects.
Fpm's user interface is modeled after Rust's Cargo, so if you're familiar with that tool,
you will feel at home with fpm.
Fpm's long term vision is to nurture and grow the ecosystem of modern Fortran applications and libraries.")
   (home-page "https://fpm.fortran-lang.org")
   (license expat)))

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
