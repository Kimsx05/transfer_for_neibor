# Runtime audit

The historical virtual-environment python symlink targets /python-env/bin/python, which is not available in this execution context. No old environment files were changed. The audit launcher uses /usr/bin/python3.12 -S with PYTHONPATH pointing to the existing cell2location environment site-packages; -S prevents unrelated system packages from leaking into this runtime. PYTHONDONTWRITEBYTECODE=1 is set for subsequent runs.

Actual package versions are in software_versions.tsv. The historical environment specification and NB model versions are copied into provenance. Current cell2location 0.1.5, scvi-tools 1.3.3 and torch 2.6.0+cu126 match the historical core versions; other package versions are recorded as actually installed, not asserted historically identical. RegressionModel and Cell2location imports and public API introspection succeeded; signatures are in logs/runtime_API_check.log. No model was instantiated or trained.

To reproduce the API inspection: bash code/python.sh code/08_runtime_api.py (from this output directory). GPU was idle at initial and final checks. No existing jobs were terminated. Standard sandbox command launch initially failed (bwrap network setup); approved escalated execution was used for project reads and v5-only writes.
