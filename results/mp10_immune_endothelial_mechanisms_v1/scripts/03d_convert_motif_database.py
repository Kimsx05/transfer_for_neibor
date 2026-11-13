import pathlib,sys,importlib.util,json,hashlib
import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
import pyarrow.feather as pf
P=pathlib.Path('resources/mp10_regulatory_priors');src=P/'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather';dst=P/'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
assert src.exists()
spec=importlib.util.spec_from_file_location('official_converter',P/'convert_cistarget_databases_v1_to_v2.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.convert_feather_v1_to_v2_vice_versa(str(src),str(dst),compression='zstd',compression_level=6,to_version=2)
a=pf.read_table(src,memory_map=True,use_threads=False);b=pf.read_table(dst,memory_map=True,use_threads=False);assert a.num_rows==b.num_rows and a.num_columns==b.num_columns
# Verify every ranking column exactly, allowing only the official features -> motifs name change.
for n in a.column_names:
 nn='motifs' if n=='features' else n
 assert a[n].equals(b[nn]),n
h=hashlib.sha256()
with dst.open('rb') as f:
 for z in iter(lambda:f.read(1024*1024),b''):h.update(z)
(P/'motif_conversion_audit.json').write_text(json.dumps(dict(source=str(src),output=str(dst),method='official aertslab Feather format conversion',all_columns_identical=True,n_rows=a.num_rows,n_columns=a.num_columns,sha256=h.hexdigest(),bytes=dst.stat().st_size),indent=2));print('CONVERSION VERIFIED',a.num_rows,a.num_columns,flush=True)
