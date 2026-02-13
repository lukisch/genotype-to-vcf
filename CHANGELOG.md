# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2026-02-13

### Initial Public Release

- PyQt6 GUI with dark theme
- VCF v4.2 output format
- Dual reference genome support (GRCh37 / GRCh38)
- Automatic build detection via dbSNP position validation
- Automatic sex detection from Y chromosome variants
- PAR (pseudo-autosomal region) handling for correct X/Y ploidy
- NCBI dbSNP REST API integration with persistent local cache
- Optional FASTA reference download with automatic indexing
- Adaptive multi-threading (4-200 workers, targeting 70% CPU usage)
- Indel detection and handling (I/D markers)
- Internal ID (i-prefix) to rsID mapping via cache lookup
- Compatible with 23andMe and other providers using the same TSV format
