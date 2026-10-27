# Independent review findings and resolution

1. Floating-point lattice boundary exclusion: fixed before any primary inference completion. COMMOT uses round(dNN,8) and 1e-6 um numeric tolerance in all ranges; SpatialCellChat uses round(dNN,8) and native tol=1e-6. Main nominal ranges remain200/105um. Native SCC distance QA confirms all12596/66166 cross-spot links within200um, including all4112/21790 boundary links. Aborted pre-tolerance COMMOT run has no completed output; its log is preserved separately.
2. Spatial epithelial caption: common renderer accepts explicit epi_caption. All spatial circles now state epithelial components in MP10-high/low-associated anchor spots; shared TME composition. Added c2l proportion-weighted mixed-expression caveat. First SCC major and EC white previews actually opened and checked.

No source data, old scRNA result, or inferred network is altered by display masks.
