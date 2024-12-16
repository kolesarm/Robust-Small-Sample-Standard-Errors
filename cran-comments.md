## Test environments
* local Ubuntu 24.04.1 LTS install, R 4.3.3
* macbuilder macOS 13.3.1 (22E261) R4.4.0 (2024-04-24)
* win-builder, R-devel and R-release
* Github actions
  - macOS 14.7.1, R 4.4.2
  - Windows Server 2022, R 4.4.2
  - Ubuntu 22.04.5 LTS, R 4.4.2
  - Ubuntu 22.04.5 LTS, R-devel
* Rhub
  - macOS 13.7.1 R-devel (2024-12-15 r87442)
  - macOS-arm64 14.7.1, R-devel (2024-12-15 r87442)
  - ubuntu-gcc12 22.04.5 LTS R-devel (2024-12-15 r87442)
  - ubuntu-nold 22.04.5 LTS, R-devel (2024-12-15 r87442)
  - ubuntu-release 22.04.5 LTS, R 4.4.2


## R CMD check results
There were no ERRORs or WARNINGs

There was 1 NOTE:

Found the following (possibly) invalid URLs:
  URL: https://www.doi.org/10.1162/REST_a_00552
    From: NEWS.md
    Status: 403
    Message: Forbidden

This is a valid URL, I checked manually.

## Downstream dependencies
There are currently no downstream dependencies for this package
