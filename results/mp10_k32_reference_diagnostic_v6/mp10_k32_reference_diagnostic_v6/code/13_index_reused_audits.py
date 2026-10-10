from common import *
import shutil
rows=[]
for name in ['scope_exclusions.tsv.gz','patient_mapping_status.tsv','frozen_source_split_manifest.tsv','split_support_by_dataset.tsv','split_support_by_subtype.tsv','new_subtype_original_major_provenance.tsv','reused_identity_K32_retention.tsv']:
 src=S/'01_scope'/name;dst=R/'00_manifest'/('reused_'+name);shutil.copy2(src,dst);assert sha(src)==sha(dst);rows.append(dict(original_path=str(src),copied_path=str(dst),sha256=sha(dst),reuse_reason='Frozen membership/source split unchanged; original remains read-only'))
wt(pd.DataFrame(rows),R/'00_manifest/REUSED_SCOPE_SOURCES.tsv')
wt(pd.DataFrame([{'purpose':x,'seed':seed(x),'master_seed':123} for x in ['reference_fit','reference_posterior','reference_reconstruction_QC','mapping_DEV','mapping_posterior_DEV','v6_control_mapping','v6_control_posterior','v6_technical_posterior']]),R/'00_manifest/SEEDS.tsv')
