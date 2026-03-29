# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Initial changelog created.

## [v1.0.2] - 2026-03-28
### Added
- Added version pinning and security hardening to all tool Dockerfiles.
- Added explicit Cigar import workaround for cuteSV Dockerfile.
- Added ARGs for tool versions in all Dockerfiles.
- Added user and workdir security best practices to Dockerfiles.
- Added image digests for all Docker images in Dockerfiles for reproducibility and traceability. Digests are specified in each tool's Dockerfile under the `tools/` directory and centrally referenced in `nextflow.config` for workflow reproducibility.

### Changed
- Updated all tool Dockerfiles to use multi-stage builds where appropriate.
- Updated all tool Dockerfiles to use minimal base images and clean up build dependencies.
- Updated README and documentation for Docker and Nextflow usage.

### Fixed
- Fixed Cigar import issues for cuteSV on Python 3.10+.
- Fixed build failures for SyRI and cuteSV on Alpine and slim images.
- Fixed missing tags/versions for all bioinformatics tools.

## [Earlier]
- See project commit history for previous changes.
