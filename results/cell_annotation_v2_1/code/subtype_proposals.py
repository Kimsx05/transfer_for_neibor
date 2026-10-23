from common import *
import numpy as np,pandas as pd
from scipy import sparse
from expression_evidence import SUBTYPE,ANCHORS
from programs import STATE
from major_proposals import ANCHOR as MAJOR_ANCHORS

def tag(s):return s.replace('/','_')
o=pd.read_parquet(ROOT/'metadata/lineage_structure.parquet');p=pd.read_parquet(ROOT/'metadata/program_evidence.parquet');genes=pd.Index(json.loads((ROOT/'metadata/observed_marker_genes.json').read_text()));raw=sparse.load_npz(ROOT/'metadata/observed_marker_counts.npz');mapping=json.loads((ROOT/'metadata/reference.json').read_text())['major_mapping'];N=len(o)
for head in ['major','subtype']:
 q=pd.concat([pd.read_csv(ROOT/f'metadata/{sec}_SingleR_{head}.tsv.gz',sep='\t').set_index('cell_id') for sec in ['Cervilla','D1_Johnson']]).reindex(o.index);assert q.candidate.notna().all()
 for c in q:o['SingleR_'+head+'_'+c]=q[c].to_numpy()
within=pd.concat([pd.read_csv(ROOT/f'metadata/{sec}_SingleR_within_major_initial.tsv.gz',sep='\t').set_index('cell_id') for sec in ['Cervilla','D1_Johnson']]).reindex(o.index)
def counts(gs):return np.asarray(raw[:,genes.get_indexer([g for g in gs if g in genes])].sum(1)).ravel()
def detected(gs):return np.asarray((raw[:,genes.get_indexer([g for g in gs if g in genes])]>0).sum(1)).ravel()
# Subtype-specific identity anchors separate generic activation from identity.
anchors={**{k:v for k,v in SUBTYPE.items()},**ANCHORS,'Classical Mono':['FCN1','S100A12','CCR2'],'cDC1':['CLEC9A','XCR1'],'cDC2':['CD1C','FCER1A','CLEC10A'],'LAMP3 DC':['LAMP3','FSCN1'],'CD4 Helper/Conventional':['CD40LG','IL17A','IL26','IL23R','IL7R'],'CD4 Naive/Memory-like':['CCR7','TCF7','LEF1','SELL'],'iCAF':['PI16','DPT','PDGFRA','CXCL12','CFD','APOD','PLA2G2A'],'Matrix_CAF':['COL11A1','COL10A1','MMP11','POSTN','CTHRC1','LRRC15','FAP'],'rCAF':['CCL19','CCL21'],'Capillary EC':['CA4','RGCC'],'Arterial EC':['GJA5','EFNB2','HEY1'],'Lymphatic EC':['PROX1','LYVE1']}
o['proposed_subtype']='Unknown';o['subtype_evidence_source']='not_resolved';o['subtype_support']='unresolved';o['subtype_ambiguity']='none';o['within_major_SingleR_initial_candidate']='not_applicable';o['within_major_SingleR_initial_margin']=np.nan
review=[];candidate_tables=[]
for major in sorted(o.proposed_major.unique()):
 if major=='Unknown':continue
 ix=np.flatnonzero(o.proposed_major.eq(major));labels=[k for k in SUBTYPE if mapping.get(k)==major];prefix='lineage_'+tag(major)
 mark=pd.read_parquet(ROOT/f'metadata/{prefix}_markers.parquet');summary=pd.read_csv(ROOT/f'metadata/{prefix}_summary.tsv',sep='\t');joint=mark[mark.scope.eq('joint')];allclusters=summary.loc[summary.scope.eq('joint'),'cluster'].tolist()
 if major in ['Epithelial','Plasma','Mast','Erythroid','Neural_glial']:
  o.loc[o.index[ix],'subtype_evidence_source']='broad_identity_only_states_examined';continue
 st=np.column_stack([p['subtype:'+lab].to_numpy()[ix] for lab in labels]);ac=np.column_stack([counts(anchors[lab])[ix] for lab in labels]);nd=np.column_stack([detected(anchors[lab])[ix] for lab in labels]);pdect=np.column_stack([detected(SUBTYPE[lab])[ix] for lab in labels]);a=np.zeros_like(st,dtype=bool);clusterok=np.zeros_like(st,dtype=bool);sources=np.full(st.shape,'not_resolved',dtype=object)
 col=major.replace('/','_');o.loc[o.index[ix],'within_major_SingleR_initial_candidate']=within.loc[o.index[ix],col+'_candidate'].to_numpy();o.loc[o.index[ix],'within_major_SingleR_initial_margin']=within.loc[o.index[ix],col+'_margin'].to_numpy()
 sr=o.SingleR_subtype_candidate.to_numpy()[ix];unpruned=~o.SingleR_subtype_pruned.to_numpy()[ix];initial=within.loc[o.index[ix],col+'_candidate'].to_numpy();mar=within.loc[o.index[ix],col+'_margin'].to_numpy();ct=o.predicted_subtype.astype(str).to_numpy()[ix]
 for j,lab in enumerate(labels):
  # Direct program or a strong identity anchor reinforced by another measured gene or a reference match.
  direct=(nd[:,j]>=2)&(pdect[:,j]>=2)&(ac[:,j]>=3)
  anchor=(ac[:,j]>=3)&(pdect[:,j]>=2)
  ref=(ac[:,j]>=1)&(sr==lab)&unpruned&(o.n_genes.to_numpy()[ix]>30)&(o.SingleR_subtype_delta_median.to_numpy()[ix]>.06)
  # Preserve subtype identity requirements (CD8/exhaustion, Tfh vs cytotoxic activation, NK vs TCR).
  compatible=np.ones(len(ix),bool)
  if lab.startswith('CD8'):compatible=detected(['CD8A','CD8B'])[ix]>0
  if lab.startswith('CD4'):compatible=(counts(['CD4'])[ix]>0)|((initial==lab)&(mar>.025))
  if lab=='Tfh/CXCL13-like':compatible=(detected(['CD4','ICOS','PDCD1','IL21','TOX2'])[ix]>0)&(detected(['GNLY','FGFBP2'])[ix]==0)
  if lab=='NK':compatible=(detected(['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'])[ix]==0)|((counts(['GNLY','KLRD1','FGFBP2'])[ix]>2*counts(['CD3D','CD3E','TRAC'])[ix])&(initial=='NK'))
  # GC identity needs coordinated GC genes in this sparse query, not BCL6/cycling alone.
  if lab=='GC B':compatible=detected(['BCL6','AICDA','RGS13'])[ix]>=2
  for cl in allclusters:
   jj=np.flatnonzero(o.lineage_cluster_v21.to_numpy()[ix]==cl);m=joint[joint.cluster.eq(cl)].set_index('gene');ags=[g for g in anchors[lab] if g in m.index];pgs=[g for g in SUBTYPE[lab] if g in m.index]
   if not ags or not pgs:continue
   ap=m.loc[ags,'detection_fraction'];apdiff=m.loc[ags,'detection_difference'];ae=m.loc[ags,'mean_log2_difference'];npan=int((m.loc[pgs,'detection_fraction']>=.07).sum());apos=int(((ap>=.07)&(ae>.08)&(apdiff>.015)).sum());coherent=(apos>=1 and npan>=2 and len(jj)>=20 and 'unclustered' not in cl)
   # Cell-level co-occurrence check avoids declaring a program assembled from disjoint cell groups.
   cofrac=float(((nd[jj,j]>0)&(pdect[jj,j]>=2)).mean()) if len(jj) else 0
   coherent=coherent and cofrac>=.05
   if coherent:
    compatible_ref=((sr[jj]==lab)&unpruned[jj])|((initial[jj]==lab)&(mar[jj]>.04))
    clusterok[jj,j]=compatible_ref&(o.n_genes.to_numpy()[ix[jj]]>30)&((ac[jj,j]>0)|(o.major_support.to_numpy()[ix[jj]]=='strong_measured'))
   review.append(dict(lineage=major,cluster=cl,subtype=lab,n=len(jj),cluster_program_supported=coherent,anchor_positive_genes=apos,panel_detected_genes=npan,member_anchor_plus_panel_fraction=cofrac,anchor_detection=';'.join(f'{g}:{ap[g]:.3f}' for g in ags),supporting_anchor_effects=';'.join(f'{g}:{ae[g]:.3f}' for g in ags)))
  a[:,j]=(direct|anchor|ref|clusterok[:,j])&compatible;sources[:,j]=np.select([direct,anchor,ref,clusterok[:,j]],['direct_marker_program','direct_marker_anchor_combination','direct_marker_plus_SingleR','cluster_supported_expression_compatible'],default='not_resolved')
 # Compare only evidence-eligible candidates; reference is a supporting tie-break, not a veto.
 rankscore=st+np.column_stack([.2*(sr==lab)+.1*(initial==lab) for lab in labels]);rankscore[~a]=-1;rank=np.argsort(-rankscore,axis=1);best=rank[:,0];v=rankscore[np.arange(len(ix)),best];accepted=v>=0
 chosen=np.array(labels)[best];src=sources[np.arange(len(ix)),best];ambiguity=np.full(len(ix),'none',object)
 if len(labels)>1:
  second=rank[:,1];v2=rankscore[np.arange(len(ix)),second];close=(v2>=0)&((v-v2)<.25);ambiguity[close]='overlapping_subtype_programs'
  # C1Q and SPP1 coexpression is retained as a state; other nearly tied identities remain broad.
  macrophage_pair=np.isin(chosen,['C1QC TAM','SPP1 TAM','FOLR2 TAM'])&np.isin(np.array(labels)[second],['C1QC TAM','SPP1 TAM','FOLR2 TAM']);accepted[close&~macrophage_pair]=False
 idx=ix[accepted];o.loc[o.index[idx],'proposed_subtype']=chosen[accepted];o.loc[o.index[idx],'subtype_evidence_source']=src[accepted];o.loc[o.index[idx],'subtype_support']=np.where(np.char.startswith(src[accepted].astype(str),'direct_marker'),'moderate_or_strong_measured','moderate_cluster_supported');o.loc[o.index[ix],'subtype_ambiguity']=ambiguity
 print(major,pd.Series(chosen[accepted]).value_counts().to_dict(),flush=True)
# State programs are descriptive and nonexclusive; never infer malignancy or mechanisms.
state_strings=np.full(N,'none',object)
for state,gs in STATE.items():
 hit=detected(gs);cnt=counts(gs);present=(hit>=2)&(cnt>=3)
 if state=='SPP1_program':present=(counts(['SPP1'])>=1)&(counts(['GPNMB','APOC1'])>=1)
 if state.startswith('Epithelial_'):present&=o.proposed_major.eq('Epithelial').to_numpy()
 if state.startswith('Plasma_'):present&=o.proposed_major.eq('Plasma').to_numpy()
 if state.startswith('Mast_'):present&=o.proposed_major.eq('Mast').to_numpy()
 if state.startswith('Myeloid_') or state in ['C1Q_program','SPP1_program','FOLR2_program']:present&=o.proposed_major.eq('Myeloid').to_numpy()
 o['state_'+state]=present;state_strings[present]=np.where(state_strings[present]=='none',state,state_strings[present]+';'+state)
o['state']=state_strings;o['C1Q_SPP1_coexpression']=o.state_C1Q_program&o.state_SPP1_program
o['annotation_support']=np.where(o.proposed_subtype.ne('Unknown'),o.subtype_support,o.major_support);o['evidence_source']=o.major_evidence_source+';'+o.subtype_evidence_source
o['unknown_reason']=np.where(o.proposed_major.eq('Unknown'),o.major_unknown_reason,np.where(o.proposed_subtype.eq('Unknown'),np.where(o.subtype_evidence_source.eq('broad_identity_only_states_examined'),'broad_identity_supported_finer_identity_not_established',np.where(o.subtype_ambiguity.ne('none'),'competing_related_subtype_programs','major_supported_subtype_information_insufficient')),'not_applicable'))
o['quality_flag']=np.where(o.mixed_observed_programs,'mixed_observed_lineage_programs','none');o.loc[o.artifact_review_flag,'quality_flag']+=';possible_segmentation_issue';o.loc[o.low_information,'quality_flag']+=';low_information';o.loc[o.mixed_or_ambiguous,'quality_flag']+=';v2_mixed_flag_preserved'
# Secondary coherent lineage transcripts are quality evidence even when one identity dominates.
secondary=[];secondary_n=np.zeros(N,dtype='int16')
for lineage,gs in MAJOR_ANCHORS.items():
 present=(detected(gs)>=2)&(counts(gs)>=3)&o.proposed_major.ne(lineage).to_numpy()&o.proposed_major.ne('Unknown').to_numpy()
 secondary.append((lineage,present));secondary_n+=present.astype('int16')
o['secondary_lineage_signal']=secondary_n>0;o['secondary_lineage_programs']='none'
for lineage,present in secondary:
 prev=o.loc[present,'secondary_lineage_programs'];o.loc[present,'secondary_lineage_programs']=np.where(prev.eq('none'),lineage,prev+';'+lineage)
o.loc[o.secondary_lineage_signal,'quality_flag']+=';secondary_lineage_transcripts_ambient_or_boundary_possible'
o['old_major_label']=o.major_label;o['old_final_label']=o.final_label;o['global_cluster_v21']=o.expression_cluster.astype(str);o['display_label']=np.where(o.proposed_subtype.ne('Unknown'),o.proposed_subtype,np.where(o.proposed_major.eq('Unknown'),'Unknown',o.proposed_major+'_subtype_unknown'))
o['major_changed_from_v2']=o.proposed_major.ne(o.major_label.astype(str));o['subtype_changed_from_v2']=o.display_label.ne(o.final_label.astype(str));o['review_status']='computational_proposal_pending_user_review';o['review_entry']='global::'+o.global_cluster_v21+';'+o.lineage_cluster_v21
# Label dictionary and all-cell table remain compact; numeric scores are diagnostics, not correct probabilities.
o.to_parquet(ROOT/'cell_annotations.parquet',compression='zstd');pd.DataFrame(review).to_csv(ROOT/'metadata/subtype_cluster_program_review.tsv',sep='\t',index=False)
js(ROOT/'metadata/subtype_rules.json',{'anchors':anchors,'reference_labels':mapping,'note':'Alternative measured-program/strong-anchor/reference/compatible-cluster routes; all score contrasts and cutoffs are descriptive proposal heuristics, not calibrated accuracy. SingleR within-major initial scores are not separate fine-tuned runs. Near-tied related identities retained as parent, C1Q/SPP1/FOLR2 programs nonexclusive. Epithelial/plasma/mast states reviewed without inventing fine identities.','cluster_route':'same-lineage observed enrichment; anchor>=7% detection, mean-log2difference>.08 and detectiondiff>.015; 2 panel genes>=7%; actual same-member anchor+panel cooccurrence>=5%; n>=20; individual reference/marker compatibility. No spatial propagation.'})
print('ALL_CELL_PROPOSALS_COMPLETE',len(o),flush=True);print(pd.crosstab(o.proposed_major,o.section_id).to_string());print(pd.crosstab(o.proposed_subtype,o.section_id).to_string())
