from common import *
import os
ri=pd.read_csv(V/'tables/ROI_review_index.tsv',sep='\t');ri.to_csv(R/'tables/fixed_ROI_index.tsv',sep='\t',index=False)
F=R/'figures';F.mkdir(exist_ok=True);rows=[]
selection={'Cervilla':'far_UMI_hotspot','D1_Johnson':'conflict_focus','UC_SCC':'typical_information','MPBC1':'typical_information','MPBC2':'high_information'}
for s,roi in selection.items():
 src=E/s/'figures'/f'ROI_single_genes_{roi}.png';name=f'{s}_raw_evidence_{roi}.png';dest=F/name
 if src.exists():
  if not dest.exists():dest.symlink_to(os.path.relpath(src,F))
  rows.append(dict(section=s,roi=roi,figure=name,source=str(src),evidence='original RNA LogNormalize single-gene expression; inherited fixed ROI, not regenerated'))
src=R.parent/'hd_bin8_banksy_v1/D1_Johnson/figures/ROI_raw_markers_conflict_focus.png';dest=F/'D1_Johnson_raw_evidence_HE_conflict_focus.png'
if not dest.exists():dest.symlink_to(os.path.relpath(src,F))
rows.append(dict(section='D1_Johnson',roi='conflict_focus',figure=dest.name,source=str(src),evidence='existing registered H&E with raw panel detection; not independent cell identity truth'))
pd.DataFrame(rows).to_csv(R/'tables/raw_ROI_figure_sources.tsv',sep='\t',index=False)
