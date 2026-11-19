from domain import *
for s in sections:
 if not (R/s/'cells_pre_MP10.parquet').exists():run(s)
