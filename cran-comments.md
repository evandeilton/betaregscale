## Submission

Update from 2.7.4 (current CRAN version) to 3.1.1; version 3.0.0 was not
submitted to CRAN, so this update includes its changes (see `NEWS.md`,
including the breaking changes listed under 3.0.0). There are no reverse
dependencies.

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.6.1, x86_64-pc-linux-gnu
  (`R CMD check --as-cran`, with the PDF manual)

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE

  Found the following (possibly) invalid URLs:
    URL: https://acervodigital.ufpr.br/xmlui/handle/1884/86624
      Status: Error
      Message: libcurl error code 60:
        SSL certificate OpenSSL verify result: unable to get local issuer
        certificate (20)

  This is the repository page of a master's dissertation. It opens without
  SSL verification: the server does not send its intermediate certificate,
  and it may respond slowly.

* checking HTML version of manual ... NOTE

  Skipping checking HTML validation: no command 'tidy' found.

  HTML Tidy is not installed on the local machine.

## Installed size

Reported as INFO, not as a NOTE: 6.6 Mb, of which libs 5.0 Mb.
