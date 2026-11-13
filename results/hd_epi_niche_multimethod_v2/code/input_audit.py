from common import *
import h5py,hashlib,importlib.metadata as md,sys,subprocess
rows=[];genes=None
for s in SECTIONS:
 base,p=source(s)
 with h5py.File(p) as f:
  g=f['matrix'];gn=g['features/id'][:].astype(str);sy=g['features/name'][:].astype(str);bc=g['barcodes'][:].astype(str);shape=g['shape'][:]
 n=pd.read_parquet(V/s/'nodes.parquet');qc=n[n.analysis_qc];assert pd.Index(bc).is_unique and pd.Index(qc.bin_id).isin(bc).all()
 if genes is None:genes=gn
 assert np.array_equal(gn,genes)
 rows.append(dict(section=s,counts_path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,genes=len(gn),supplied_bins=len(bc),analysis_bins=len(qc),anchors=int(n.epi_anchor.sum()),common_gene_order=True,barcode_QC_coverage=1.0,original_tissue_mask_known=bool(n.in_tissue_original.notna().all()),coordinate_units='um',coordinate_direction='increasing y down; relative array only for MPBC'))
pd.DataFrame(rows).to_csv(R/'tables/counts_input_audit.tsv',sep='\t',index=False)
versions={p:md.version(p) for p in ['numpy','scipy','pandas','anndata','scanpy','squidpy','scikit-learn','igraph','leidenalg','harmonypy','MENDER','cellcharter','scvi-tools','torch','torchgmm']};versions['python']=sys.version;versions['MENDER_commit']=subprocess.check_output(['git','-C',str(R/'environment/MENDER'),'rev-parse','HEAD'],text=True).strip();(R/'versions.json').write_text(json.dumps(versions,indent=2));print(pd.DataFrame(rows)[['section','analysis_bins','anchors','genes']].to_string(index=False))
