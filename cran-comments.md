## Submission

Update of `betaregscale` to version 3.1.0. The changes are listed in
`NEWS.md`.

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.6.1, x86_64-pc-linux-gnu
  (`R CMD check --as-cran`, with the PDF manual)

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE

  Found the following (possibly) invalid URLs:
    URL: https://hdl.handle.net/1884/86624
      Status: Error
      Message: libcurl error code 60:
        SSL certificate OpenSSL verify result: unable to get local issuer
        certificate (20) (Status without verification: OK)

  This is the persistent handle of a master's dissertation. It redirects to
  `acervodigital.ufpr.br`, whose server does not send its intermediate
  certificate; the link opens in web browsers.

* checking HTML version of manual ... NOTE

  Skipping checking HTML validation: no command 'tidy' found.

  HTML Tidy is not installed on the local machine.

## Installed size

Reported as INFO, not as a NOTE: 7.5 Mb, of which libs 5.0 Mb. It can be
larger on Linux builds that add `-g` to `CXXFLAGS`, since the shared library
then carries debug information.

## Reverse dependencies

There are no reverse dependencies.
