## Submission

Update of `betaregscale`. The changes are listed in `NEWS.md`.

## Test environments

* Local: Ubuntu (Linux 7.0.0), R 4.6.1, x86_64-pc-linux-gnu
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1),
  macOS-latest (R release), windows-latest (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* The installed size can be reported as large on Linux builds that add `-g`
  to `CXXFLAGS`: the shared library then carries debug information (about
  325 KB when stripped).
* `CXX_STD = CXX17` is set because the package supports R >= 4.1.0, which
  does not default to C++17.

## Reverse dependencies

There are no reverse dependencies.
