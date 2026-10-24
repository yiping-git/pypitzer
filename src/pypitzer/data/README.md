# Runtime thermodynamic database

`parameters-v2026.1.sqlite` is the curated runtime snapshot extracted from the
private reference implementation for PyPitzer 2.0.0. It contains 96 binary
records, 123 ternary records, and 70 reaction records. Runtime code opens it in
SQLite `mode=ro` with `PRAGMA query_only = ON`.

- Database version: `2026.1`
- SHA-256: `13a5c8ac3d6323218898ca85758f3ff67ba3235cd90148c492e4abf7aae2eadf`

No database creation, editing, literature-import, or curation tooling is
distributed with this package.
