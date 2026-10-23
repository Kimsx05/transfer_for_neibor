from common import *
import numpy as np,pandas as pd
from scipy import sparse
from programs import MAJOR
# Identity anchors avoid using ubiquitous ECM/contractile/secretory genes as lineage proof.
ANCHOR={'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7','KRT5','KRT14','UPK1A','UPK2'],'T/NK':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','GNLY','KLRD1'],'B':['CD79A','CD79B','MS4A1','CD19','CD22'],'Plasma':['MZB1','JCHAIN','DERL3','TNFRSF17'],'Myeloid':['LST1','TYROBP','FCER1G','CSF1R','CD68','AIF1','C1QA','C1QB','C1QC'],'Fibroblast':['DCN','LUM','PDGFRA','DPT'],'Mural':['RGS5','CSPG4','PDGFRB','MYH11','CNN1','MYOCD','LMOD1'],'Endothelial':['PECAM1','VWF','EMCN','KDR','CLDN5','ESAM','PROX1'],'Mast':['TPSAB1','TPSB2','CPA3','KIT','MS4A2'],'Erythroid':['HBB','HBA1','HBA2','ALAS2','GYPA','KLF1'],'Neural_glial':['SOX10','S100B','PLP1','MPZ','GFAP']}
def run():
 o=pd.read_parquet(ROOT/'metadata/cell_metadata.parquet');p=pd.read_parquet(ROOT/'metadata/program_evidence.parquet');genes=pd.Index(json.loads((ROOT/'metadata/observed_marker_genes.json').read_text()));raw=sparse.load_npz(ROOT/'metadata/observed_marker_counts.npz');names=list(MAJOR);N=len(o);ac=[];ad=[]
 for label in names:
  r=raw[:,genes.get_indexer([g for g in ANCHOR[label] if g in genes])];ac.append(np.asarray(r.sum(1)).ravel());ad.append(np.asarray((r>0).sum(1)).ravel())
 ac=np.array(ac).T;ad=np.array(ad).T;v=p[['major:'+m for m in names]].to_numpy(copy=True);v[ac==0]=0
 # Same-anchor weighted log-expression ranks alternatives. Broad matrix alone cannot establish fibroblast.
 rank=np.argsort(-v,axis=1);best=rank[:,0];second=rank[:,1];i=np.arange(N);s=v[i,best];s2=v[i,second];dom=s/(s2+.35);hits=ad[i,best];counts=ac[i,best];cand=np.array(names)[best]
 sr=[]
 for sec in ['Cervilla','D1_Johnson']:
  q=pd.read_csv(ROOT/f'metadata/{sec}_SingleR_major.tsv.gz',sep='\t').set_index('cell_id');sr.append(q)
 sr=pd.concat(sr).reindex(o.index);assert sr.candidate.notna().all()
 # Separate routes: direct coherent measured signal; a strong anchor; expression-compatible reference; or measured signal supported by a global expression group.
 strong=((hits>=3)&(dom>=1.35)&(s>=1.2))|((hits>=2)&(counts>=2)&(dom>=1.1)&(ad[i,second]<2)&(o.n_genes.to_numpy()>30))
 anchor=(counts>=4)&(dom>=1.8)&(o.n_genes.to_numpy()>30)
 reference=(counts>=1)&(dom>=1.5)&(sr.candidate.to_numpy()==cand)&(~sr.pruned.to_numpy())&(sr.delta_median.to_numpy()>.05)&(o.n_genes.to_numpy()>30)
 groupok=np.zeros(N,bool);globalcontext={}
 for cluster,idx in o.groupby('expression_cluster',observed=True).indices.items():
  idx=np.array(idx);direct=idx[strong[idx]|anchor[idx]];freq=pd.Series(cand[direct]).value_counts(normalize=True)
  globalcontext[str(cluster)]={'n':len(idx),'direct_expression_composition':freq.to_dict()}
  if str(cluster)=='not_in_structure_core' or not len(freq):continue
  lead=freq.index[0]
  # New measured-expression composition, never old candidate purity. This is only one supporting route.
  if freq.iloc[0]>=.7:
   groupok[idx]=(cand[idx]==lead)&(counts[idx]>=1)&(dom[idx]>=1.8)&(o.n_genes.to_numpy()[idx]>30)&((sr.candidate.to_numpy()[idx]==lead)|(o.predicted_major.astype(str).to_numpy()[idx]==lead))
 assigned=strong|anchor|reference|groupok
 # Sparse cells can still have a strong multi-gene observation; no low-info deletion.
 out=np.where(assigned,cand,'Unknown').astype(object);source=np.select([strong,anchor,reference,groupok],['direct_marker_program','direct_marker_anchor','direct_marker_plus_SingleR','cluster_supported_compatible_marker'],default='insufficient_evidence').astype(object)
 # Comparable coherent incompatible programs within one object: flag, do not call a new type/doublet.
 comparable=(ad[i,best]>=2)&(ad[i,second]>=2)&(ac[i,best]>=3)&(ac[i,second]>=3)&(s2>=1.2)&(dom<1.45)
 # T/B and B/plasma transitional possibilities are unresolved, not presumed segmentation errors.
 mixed=comparable;out[mixed]='Unknown';source[mixed]='conflicting_observed_programs'
 # Reference-external erythroid/glial require coordinated non-generic evidence and normal depth.
 ext=np.isin(cand,['Erythroid','Neural_glial']);extok=ext&(hits>=3)&(counts>=6)&(dom>=1.8)&(o.n_genes.to_numpy()>30);out[ext&~extok]='Unknown';source[ext&~extok]='reference_external_signal_insufficient'
 morph=o.morphology_review_flag.astype(str).str.contains('ambiguous|possible|irregular',case=False,regex=True).to_numpy();out[morph]='Unknown';source[morph]='existing_morphology_ambiguity'
 # Compatible sparse members can receive expression-projection support from independently marker-defined seeds.
 from pynndescent import NNDescent
 coords=np.load(ROOT/'metadata/existing_coordinates.npz');latent=coords['X_scVI'];valid=np.isfinite(latent).all(1);seed=valid&(out!='Unknown')&(o.n_genes.to_numpy()>30)&~ext
 seedix=np.flatnonzero(seed);queryix=np.flatnonzero(valid&(out=='Unknown')&~mixed&~morph&(o.n_genes.to_numpy()>30)&~ext)
 graph=NNDescent(latent[seedix],n_neighbors=31,random_state=123,n_jobs=8);graph.prepare();ni,dd=graph.query(latent[queryix],k=30);sl=out[seedix];lab=sl[ni];proposal=[];pur=[]
 for ll in lab:
  uu,cc=np.unique(ll,return_counts=True);kk=cc.argmax();proposal.append(uu[kk]);pur.append(cc[kk]/len(ll))
 proposal=np.array(proposal);pur=np.array(pur);pi=np.array([names.index(k) for k in proposal]);own=ac[queryix,pi];ownn=ad[queryix,pi];other=ac[queryix].copy();other[np.arange(len(queryix)),pi]=0;othern=ad[queryix].copy();othern[np.arange(len(queryix)),pi]=0
 radii=graph.neighbor_graph[1][:,-1];radius_by={lab:float(np.quantile(radii[sl==lab],.95)) for lab in np.unique(sl)};close=dd[:,-1]<=np.array([radius_by[k] for k in proposal])*1.2
 compatible=(othern.max(1)<2)|(own>=2*other.max(1));reference_match=(sr.candidate.to_numpy()[queryix]==proposal)|(o.predicted_major.astype(str).to_numpy()[queryix]==proposal)
 project=(pur>=.85)&close&compatible&((own>0)|((o.n_genes.to_numpy()[queryix]>=60)&reference_match))
 qi=queryix[project];out[qi]=proposal[project];source[qi]='expression_projection_compatible';projection_purity=np.full(N,np.nan);projection_purity[queryix]=pur;o['expression_projection_neighbor_fraction']=projection_purity
 js(ROOT/'metadata/major_projection.json',dict(seed_cells=len(seedix),query_cells=len(queryix),assigned=len(qi),method='Existing integrated scVI20; NNDescent30 expression neighbors of independently marker-supported seeds; seed123; >=85% new seed-label consistency, within 1.2x lineage seed95% radius, no coherent competing anchors. No physical neighbor propagation.',radius_by_lineage=radius_by))
 print('PROJECTION',len(seedix),len(queryix),len(qi),flush=True)
 o['proposed_major']=out;o['major_evidence_source']=source;o['major_support']=np.select([source=='expression_projection_compatible',strong&assigned,anchor&assigned,reference|groupok],['moderate_expression_projection','strong_measured','moderate_measured','moderate_supported'],default='unresolved');o.loc[o.proposed_major.eq('Unknown'),'major_support']='unresolved'
 o['major_top_program']=cand;o['major_competing_program']=np.array(names)[second];o['major_program_ratio']=dom;o['major_anchor_counts']=counts;o['major_anchor_genes']=hits;o['mixed_observed_programs']=mixed;o['artifact_review_flag']=morph|((cand=='Erythroid')&(counts>=4));o['artifact_evidence']=np.select([morph,(cand=='Erythroid')&(counts>=4)],['existing_v2_morphology_ambiguity','hemoglobin_dominant_object_check_nuclear_segmentation'],default='none')
 o['major_unknown_reason']=np.select([out!='Unknown',morph,mixed,ext&~extok,o.low_information.to_numpy(),counts==0],['not_applicable','possible_segmentation_issue','competing_lineage_programs','reference_external_signal_insufficient','low_information_identity_unresolved','no_identity_anchor_detected'],default='sparse_or_competing_expression')
 o.to_parquet(ROOT/'metadata/major_proposals.parquet',compression='zstd');js(ROOT/'metadata/global_measured_composition.json',globalcontext);js(ROOT/'metadata/major_rules.json',dict(anchors=ANCHOR,note='Heuristic review proposals, not calibrated probabilities. Alternative direct-marker, strong-anchor, reference-supported and measured-cluster-compatible routes. Old calibration and tool agreement are never mandatory. Ratios are program contrasts; missing markers alone are not contradictory expression.',strong='3 observed anchors,ratio>=1.35,score>=1.2 OR 2 specific anchors with no coherent competitor,ratio>=1.1,n_genes>30',anchor='>=4 anchor UMIs,ratio>=1.8,n_genes>30',reference='any anchor,ratio>=1.5,SingleR unpruned delta>.05,n_genes>30',cluster='measured coherent global context (>=70% direct-expression members), any anchor,ratio>=1.8,one compatible classifier,n_genes>30',mixed='two coherent competing programs of comparable weighted strength',external='3 anchors,6UMIs,ratio>=1.8,n_genes>30; manual review remains required'))
 print(pd.crosstab(o.proposed_major,o.section_id).to_string(),flush=True)
if __name__=='__main__':run()
