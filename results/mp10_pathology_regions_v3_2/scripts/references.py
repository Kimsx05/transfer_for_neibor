from run_v3_2 import *
# Published S1 panels only. Cropping is for display, never segmentation or barcode label transfer.
layout=[['BC01','BC11','BC02'],['BC05','BC06','BC07'],['BC08','BC09','BC10'],['BC12','BC03','BC04'],['BC15','BC13','BC16'],['BC17','BC18','BC14'],['BC19','BC20','BC22'],['BC21']]
s1=Image.open(O/'sources/image1.png');positions={s:(r,c) for r,row in enumerate(layout) for c,s in enumerate(row)}
notes={
'BC01':'S1 contains Tumor Cell/Tumor Stroma in the left/central tissue and extensive Normal Stroma on the right/lower side. Fragmentation and holes are visible. Molecular-only boundary cannot establish the exact macro transition; unresolved areas are retained.',
'BC18':'S1 central/lower blue Tumor Stroma and salmon Tumor Cell overlap qualitatively with former outside. D8 contextual extension is calibrated to coexistence, not pixel-transferred truth. Gray Others remain unassigned; no whole-section-inside claim.',
'BC13':'S1 has substantial Normal Stroma, Muscle and Necrosis plus Tumor Cell/Tumor Stroma. A single high-epithelial mask is incomplete; exact normal/tumour boundary cannot be transferred from figure.',
'BC11':'S1 shows interspersed Tumor Cell and Tumor Stroma, with gray Others. Low malignant composition is not necessarily macro outside.',
'BC02':'S1 shows Tumor Cell, Immune Cell and gray Others across separated pieces. Immune-rich tissue is not automatically outside.',
'BC05':'S1 shows a tumour cell/stroma focus with extensive Others. Others is neither verified normal nor lesion.',
'BC06':'S1 shows tumour cells interspersed with muscle and stroma; muscle alone is not independent proof of an invasive front.',
'BC07':'S1 separates tumour-associated tissue from a lower normal-stromal compartment; low-resolution image cannot define barcode assignments.',
'BC08':'S1 includes broad muscle and necrotic compartments alongside tumour cells; no automatic inclusion/exclusion of all those compartments.',
'BC09':'S1 shows a small tumour-cell-containing tissue fragment; capture/tissue cut is not tumour outer boundary.',
'BC10':'S1 shows a tumour-cell-dominated round tissue section and a central void; neither the round cut edge nor the void is automatically a macro lesion edge.',
'BC12':'S1 includes a large Tumor Stroma compartment adjoining Tumor Cell; treating stroma as outside would narrow lesion extent.',
'BC03':'S1 includes tumour cells and necrosis with surrounding Others; exact extent remains qualitative.',
'BC04':'S1 has multiple separated tissue pieces with tumour-associated stroma/cells and muscle; no cross-fragment propagation.',
'BC15':'S1 includes extensive necrosis and Tumor Stroma; molecular low-signal areas need independent extent treatment.',
'BC16':'S1 shows separated, often Tumor Stroma-rich fragments; low malignant abundance alone cannot classify entire pieces as outside.',
'BC17':'S1 shows elongated separate tumour-containing pieces with stroma/muscle; no blanket convex-hull joining.',
'BC14':'S1 shows extensive Tumor Stroma and smaller tumour-cell/necrotic regions. S9 assigns broad blue normal-distance classes, illustrating that S1 compartments and S9 zonal labels are not interchangeable macro truth.',
'BC19':'S1 contains tumour foci amid Others; missing annotation is not negative lesion evidence.',
'BC20':'S1 shows scattered tumour cells and immune aggregates among Others; no exact mask inferred from the picture.',
'BC21':'S1 shows interspersed tumour cells, stroma and immune aggregates; composition boundaries need not be macro outer boundaries.',
'BC22':'S1 shows tumour cell/stroma and immune foci amid Others; exact macro extent unavailable.'}
rows=[]
for sec in M.index:
 dataset=M.loc[sec,'dataset_id'];out=O/'by_dataset'/dataset/sec
 row={'dataset_id':dataset,'section_id':sec,'reference_form':'unavailable','annotation_level':'none','sample_mapping':'not established','used_for_rule_calibration':False,'quantitative_overlap':'not computed; no matched barcode labels','reference_url':'','comparison':''}
 if sec.startswith('GSE319536_'):
  sample=sec.split('_')[-1];row.update(reference_form='published S1/S9 raster figures + author R script; no public per-spot RDS in repository tree',annotation_level='S1 Level2 compartments; author code Level1 Tumor/Normal/Other Region; S9 relative distance zones',sample_mapping='exact displayed BC alias; qualitative same-section comparison',reference_url='https://doi.org/10.1158/2159-8290.CD-26-0099',comparison=notes.get(sample,'Same BC alias in S1; qualitative comparison only.'),used_for_rule_calibration=sec in CFG['author_tumour_stroma_context'])
  if sample in positions and (out/'complete.json').exists():
   rr,cc=positions[sample];crop=s1.crop((80+350*cc,25+202*rr,410+350*cc,min(1658,221+202*rr)))
   d=pd.read_parquet(out/'spot_regions_v3_2.parquet');fig,axs=plt.subplots(1,3,figsize=(15,6))
   for ax,c,title in zip(axs[:2],['v3_1_lesion_location','lesion_location'],['v3.1 (project coordinates)','v3.2 (same project coordinates)']):
    ax.scatter(d.x,d.y,s=4,c=[LP.get(v,'#999999') for v in d[c]],marker='h',linewidths=0);ax.invert_yaxis();ax.set_aspect('equal');ax.axis('off');ax.set_title(title)
   axs[2].imshow(crop);axs[2].axis('off');axs[2].set_title('Published S1 annotation + H&E\nSeparate figure scale; NOT registered barcode labels')
   legend_ax=fig.add_axes([.68,.10,.29,.10]);legend_ax.imshow(s1.crop((450,1490,1130,1590)));legend_ax.axis('off')
   axs[0].legend(handles=[Patch(color='#bd3345',label='Project inside'),Patch(color='#4389b6',label='Project outside')],fontsize=7,loc='lower center',bbox_to_anchor=(.5,-.12),frameon=False);axs[1].legend(handles=[Patch(color=LP[k],label=LL[k]) for k in LP if k in set(d.lesion_location)],fontsize=7,loc='lower center',bbox_to_anchor=(.5,-.12),frameon=False)
   fig.suptitle(sec+' | author comparison is qualitative, not segmentation accuracy');fig.text(.05,.025,row['comparison'],fontsize=7,wrap=True);fig.savefig(out/'05_author_S1_comparison.png',dpi=300,bbox_inches='tight');fig.savefig(out/'05_author_S1_comparison.pdf',bbox_inches='tight');plt.close(fig)
 elif sec.startswith('GSE285715_'):
  row.update(reference_form='published marker-based cell-region figures + GEO sample metadata',annotation_level='tumour/fibroblast/immune cell composition; not macro lesion outer truth',sample_mapping='local T1-T3/N1-N3 GEO aliases confirmed; no verified Tumor-419a to T3 crosswalk',reference_url='https://doi.org/10.1073/pnas.2505125123',comparison='Do not equate fibroblast/immune regions with outside. T3 widespread EPI and tumour source do not establish whole-field inside; no author barcode macro label available.')
 elif sec.startswith('NPJ_'):
  row.update(reference_form='paper tumour/interface figures; local alias bridge unavailable',annotation_level='SpaCET tumour and invasive-front analysis',reference_url='https://doi.org/10.1038/s41698-026-01417-x',comparison='sample1/sample2 mapping to paper/GEA not verified. Shared SpaCET/expression evidence is not independent validation; no transferred labels or invasive-front claim.')
 else:
  row['comparison']='No usable same-section professional macro-region annotation obtained in existing project sources; retained unverified identities. No validation pass inferred from missing reference.'
 rows.append(row)
pd.DataFrame(rows).to_csv(O/'tables/reference_comparison.tsv',sep='\t',index=False)
