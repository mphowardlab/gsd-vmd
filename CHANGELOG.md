# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2019-02-13
### Changed
- The embedded GSD library has been updated to v1.6.0.
- The license year has been updated.
- gsd-vmd is now hosted by Michael Howard on github instead of bitbucket.

## [0.1.2] - 2017-09-26
### Fixed
- Fixed a CMake cache error when VMD is not found.

## [0.1.1] - 2017-01-16
### Fixed
- Ensured unit tests are compatible with compilers without C++11 support.
- Fixed rpath for testing so that directory structure is not assumed.

## 0.1.0 - 2017-01-12
### Added
- Implementation of the VMD molfile plugin for HOOMD-blue GSD files.
The plugin is able to read particle and bond data. Angles, dihedrals,
and writing are not currently supported.

[Unreleased]: https://github.com/mphoward/gsd-vmd/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/mphoward/gsd-vmd/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/mphoward/gsd-vmd/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/mphoward/gsd-vmd/compare/v0.1.0...v0.1.1
