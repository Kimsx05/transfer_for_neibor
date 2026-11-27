from base_data import *
import datetime
check=pd.read_csv(O/'tables/delivery_checks.tsv',sep='\t');assert check['pass'].all()
rows=[
 ('A1','completed','Cervilla,D1_Johnson','Original 2um counts, verified primary and image transform, full retained-primary raw count comparison.'),
 ('A1-other','missing_input','UC_SCC,MPBC1,MPBC2','Complete local 2um GEX+positions+registered H&E not located.'),
 ('A2','completed','HD two sections','Ten fixed buffered ROIs, old and exact HE-only reuse vs new secondary. Morphology examples examined.'),
 ('A3','completed','HD two sections','Frozen H&E-priority GEX segmentation; secondary optical QC updated after whole-section examples.'),
 ('A4','completed','HD two sections','Same-cell CellTypist vs TACCO OT fractions, lineage cluster expression and original-ID decisions.'),
 ('A5','completed','HD two sections','Raw-count reconstruction, common scVI, within-section graph and one shared CellCharter.'),
 ('B1','completed','Xenium five sections','Original segmentation, panel, old scores, individual original gate reasons.'),
 ('B2','completed','Xenium five sections','Specimen/patient-group held-out results, panel-restricted BLCA hierarchical CellTypist.'),
 ('B3','completed','Xenium five sections','One shared scVI using original cells and counts.'),
 ('B4','completed','Xenium five sections','Lineage clusters, four levels, marker dotplots and spatial/segmentation review.'),
 ('B5','completed','Xenium five sections','One shared CellCharter model, balanced training after full graph aggregation.')]
write(pd.DataFrame(rows,columns=['module','status','scope','evidence']),'module_status.tsv')
scope=Path('/sys/fs/cgroup/user.slice/user-1291.slice/user@1291.service/blca.slice/blca-retry27.slice')
stats={k:(scope/k).read_text().strip() for k in ['memory.max','memory.peak','pids.max','pids.peak','memory.events'] if (scope/k).exists()}
main=[]
for line in (O/'logs/resources.tsv').read_text().splitlines():
 v=line.split('\t')
 if len(v)==4:main.append((v[0].split('code/')[-1],float(v[1]),int(v[2]),int(v[3])))
a=pd.DataFrame(main,columns=['stage','seconds','peak_RSS_KiB','return_code']);write(a,'stage_resources.tsv')
(O/'logs/resource_summary.json').write_text(json.dumps({'shared_slice_after_setup':stats,'maximum_instrumented_stage_RSS_KiB':int(a.peak_RSS_KiB.max()),'logical_CPU_affinity':[2,3,4,5],'device_used':'CPU','all_observed_initial_scopes':'7-17 tasks per heavy scope, same four permitted compute CPUs. Initial cgroups separate; no continuous aggregate peak before shared slice.','scope_error_handling':'initial TensorFlow attempt abort signal 6 and initial HD scVI exit 143 separately logged. Successful replacements listed in runtime_notes.'},indent=2))

