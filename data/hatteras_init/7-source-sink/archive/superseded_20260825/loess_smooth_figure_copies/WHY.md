# Superseded figure copies

These sat in the calibration output folder and were OLDER than the
figures in `7-source-sink/figures/`, which is where the analysis scripts
actually write. `export_be_calibration.py` copied these over the current
ones every time it ran, so the record silently reverted.

Moved here 2026-09-12. `figures/` is the record now, and the exporter
reads the same folder rather than copying into it.
