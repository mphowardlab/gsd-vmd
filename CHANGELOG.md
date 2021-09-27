# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.1] - 2021-09-27
### Fixed
- Set triclinic box size correctly (bug fix thanks to Kevin Silmore).

### Changed
- Update code of conduct. The default git branch is renamed `main`. More
  [information](https://sfconservancy.org/news/2020/jun/23/gitbranchname) is available.

## [0.4.0] - 2021-04-12
### Changed
- The embedded GSD library has been updated to v2.4.1.
- The license has been updated. gsd-vmd is now maintained as part of our work
  at Auburn University.

## [0.3.0] - 2020-02-13
### Changed
- The embedded GSD library has been updated to v2.0.0. This version supports
reading both GSD 1.0 (HOOMD < 2.9) and GSD 2.0 (HOOMD >= 2.9) files.
- CMake can optionally configure 64-bit-only builds for macOS.
- The license year has been updated.

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

[Unreleased]: https://github.com/mphowardlab/gsd-vmd/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/mphowardlab/gsd-vmd/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/mphowardlab/gsd-vmd/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/mphowardlab/gsd-vmd/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/mphowardlab/gsd-vmd/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/mphowardlab/gsd-vmd/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/mphowardlab/gsd-vmd/compare/v0.1.0...v0.1.1
