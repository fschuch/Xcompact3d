# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

1. CHANGELOG.md by [@fschuch](https://github.com/fschuch).
1. Xdmf writer included back, making easy to visualize the binary fields in other tools, by [@fschuch](https://github.com/fschuch).
1. Sandbox flow configuration ([BC-Sandbox.f90](./src/BC-Sandbox.f90)), by [@fschuch](https://github.com/fschuch).
1. Immersed boundary method for Scalar field(s), together with a new Lagrangian Interpolator for no-flux BC at solid/fluid interface, by [@fschuch](https://github.com/fschuch).
1. Module Probes ([tools.f90](./src/tools.f90)), so now it can be used for any flow configuration, by [@fschuch](https://github.com/fschuch).

### Changed

1. Enumeration, directory structure and filename format for binary outputs by [@fschuch](https://github.com/fschuch), details in [#3](https://github.com/fschuch/Xcompact3d/issues/3).
1. Subroutine `read_geomcomplex` [genepsi3d.f90](./src/genepsi3d.f90) was modified ir order to be compatible with the new [xcompact3d_toolbox](https://github.com/fschuch/xcompact3d_toolbox), by [@fschuch](https://github.com/fschuch).
1. README was recovered from previous public repository, by [@fschuch](https://github.com/fschuch).
1. CI changed from Travis to GitHub actions, by [@fschuch](https://github.com/fschuch).

### Removed

1. `subroutine VISU_PRE` ([visu.f90](./src/visu.f90)) is obsolete, since pressure is written at `subroutine write_snapshot` ([visu.f90](./src/visu.f90)), by [@fschuch](https://github.com/fschuch).

### Fixed

1. Calling `write_snapshot` and `postprocess_case` for `itime=0`, by [@fschuch](https://github.com/fschuch).
1. [#4](https://github.com/fschuch/Xcompact3d/issues/4) Wrong derivative routine (and Boundary Conditions) for scalar field, by [@fschuch](https://github.com/fschuch).
1. [#6](https://github.com/fschuch/Xcompact3d/issues/6) Not computing gravitational term, by [@fschuch](https://github.com/fschuch).

[Unreleased]: https://github.com/xcompact3d/Incompact3d/compare/dev...fschuch:master
