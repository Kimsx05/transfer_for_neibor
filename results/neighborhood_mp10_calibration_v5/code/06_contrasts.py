from common import *
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
SCORES=['spot_MP10_model_universe','epi_MP10_model_universe','spot_MP10_frozen']

def main(out):
 P=out/'04_fixed_niche_contrasts';G=out/'02_gene_evidence';allm=pd.read_parquet(G/'spatial_technical_metadata.parquet');d=allm[allm.accepted_label.ge(0)&allm.epi_score_valid].copy();man=pd.read_csv(V4/'03_epi_attribution/model_and_gene_universe_manifest.tsv',sep='\t');rows=[];match=[];coverage=[];gene=[];loo=[];rng=np.random.default_rng(123)
 for mr in man.itertuples():
  g=d[d.section_id.eq(mr.section_id)].copy();ar=np.load(G/'by_section'/f'{mr.section_id}.npz') if (G/'by_section'/f'{mr.section_id}.npz').exists() else None;axis=pd.Index(ar['spot_uid']) if ar is not None else None
  for A,B in PAIRS:
   name=f'niche{A}_minus_niche{B}';aa=g[g.accepted_label.eq(A)];bb=g[g.accepted_label.eq(B)];valid=len(aa)>=20 and len(bb)>=20;coverage.append(dict(section_id=mr.section_id,dataset_id=mr.dataset_id,pair=name,n_A=len(aa),n_B=len(bb),status='ELIGIBLE' if valid else 'INSUFFICIENT_SAME_SECTION_FULL20',reason='' if valid else mr.status))
   if not valid:continue
   base=dict(section_id=mr.section_id,dataset_id=mr.dataset_id,analysis_unit_uid=g.analysis_unit_uid.iloc[0],analysis_unit_type=g.analysis_unit_type.iloc[0],known_patient=bool(g.source_patient_identity_known.iloc[0]),pair=name)
   sets={'MAIN':(aa,bb),'EPI_GE_0.5':(aa[aa.epi_fraction_center_mean.ge(.5)],bb[bb.epi_fraction_center_mean.ge(.5)]),'EPI_GE_0.7':(aa[aa.epi_fraction_center_mean.ge(.7)],bb[bb.epi_fraction_center_mean.ge(.7)]),'COMPLETE_RINGS':(aa[aa.complete_rings],bb[bb.complete_rings])}
   # Greedy nearest technical match without replacement, no outcome in matching.
   sd=np.std(np.log1p(pd.concat([aa,bb]).model_UMI),ddof=1);cal=.25*sd;ma=[];mb=[]
   if cal>0:
    Xb=np.column_stack([bb.epi_fraction_center_mean/.05,np.log1p(bb.model_UMI)/cal]);Xa=np.column_stack([aa.epi_fraction_center_mean/.05,np.log1p(aa.model_UMI)/cal]);tree=cKDTree(Xb);used=set()
    for i in rng.permutation(len(aa)):
     cand=[j for j in tree.query_ball_point(Xa[i],r=1,p=np.inf) if j not in used]
     if cand:
      j=min(cand,key=lambda j:(float(np.sum((Xb[j]-Xa[i])**2)),j));used.add(j);ma.append(i);mb.append(j)
   am=aa.iloc[ma];bm=bb.iloc[mb];sets['TECH_MATCHED']=(am,bm);match.append(dict(**base,n_before_A=len(aa),n_before_B=len(bb),n_pairs=len(ma),epi_caliper=.05,logUMI_caliper=cal,epi_mean_diff_before=aa.epi_fraction_center_mean.mean()-bb.epi_fraction_center_mean.mean(),epi_mean_diff_after=am.epi_fraction_center_mean.mean()-bm.epi_fraction_center_mean.mean(),logUMI_mean_diff_before=np.log1p(aa.model_UMI).mean()-np.log1p(bb.model_UMI).mean(),logUMI_mean_diff_after=np.log1p(am.model_UMI).mean()-np.log1p(bm.model_UMI).mean(),status='CALCULABLE' if len(ma)>=20 else 'MATCHED_COVERAGE_INSUFFICIENT'))
   for aset,(a,b) in sets.items():
    ok=len(a)>=20 and len(b)>=20
    for score in SCORES:rows.append(dict(**base,analysis_set=aset,score=score,n_A=len(a),n_B=len(b),effect=a[score].median()-b[score].median() if ok else np.nan,status='CALCULABLE' if ok else 'INSUFFICIENT_AFTER_FILTER',epi_A=a.epi_fraction_center_mean.median(),epi_B=b.epi_fraction_center_mean.median(),UMI_A=a.model_UMI.median(),UMI_B=b.model_UMI.median(),nFeature_A=a.model_nFeature.median(),nFeature_B=b.model_nFeature.median(),complete_A=a.complete_rings.mean(),complete_B=b.complete_rings.mean(),pathology_scope='unknown'))
   for score in SCORES:
    av=aa.groupby([aa.array_row//8,aa.array_col//8])[score].median();bv=bb.groupby([bb.array_row//8,bb.array_col//8])[score].median();rows.append(dict(**base,analysis_set='EQUAL_SPATIAL_BLOCK_MEDIANS',score=score,n_A=len(aa),n_B=len(bb),n_blocks_A=len(av),n_blocks_B=len(bv),effect=av.median()-bv.median(),status='DESCRIPTIVE_BLOCK_SENSITIVITY_NOT_PATIENTS'))
   ia=axis.get_indexer(aa.spot_uid);ib=axis.get_indexer(bb.spot_uid);den=float(g.AUC_denominator.iloc[0]);K=int(np.ceil(.05*g.n_model_genes.iloc[0]));den19=sum(K-j for j in range(1,20))
   for layer,col in [('obs','spot_MP10_model_universe'),('epi','epi_MP10_model_universe')]:
    for j,gg in enumerate(TOP):
     c=ar[layer+'_contrib'][:,j];gene.append(dict(**base,layer=layer,gene=gg,auc_mean_difference=c[ia].mean()-c[ib].mean(),counts_median_difference=np.median(ar[layer+'_counts'][ia,j])-np.median(ar[layer+'_counts'][ib,j]),rank_median_difference=np.median(ar[layer+'_rank'][ia,j])-np.median(ar[layer+'_rank'][ib,j])))
     sa=aa[col].to_numpy()-c[ia];sb=bb[col].to_numpy()-c[ib];loo.append(dict(**base,layer=layer,removed_gene=gg,original_effect=aa[col].median()-bb[col].median(),fixed_denominator_effect=np.median(sa)-np.median(sb),renormalized19_effect=(np.median(sa)-np.median(sb))*den/den19,n_A=len(ia),n_B=len(ib)))
 write_tsv(out,P/'pair_coverage_all_sections.tsv',coverage);ef=pd.DataFrame(rows);write_tsv(out,P/'paired_section_effects.tsv',ef);write_tsv(out,P/'matching_balance_and_effects.tsv',match)
 u=ef.dropna(subset=['effect']).groupby(['analysis_set','score','pair','dataset_id','analysis_unit_uid','analysis_unit_type','known_patient'],dropna=False).agg(effect=('effect','mean'),n_sections=('section_id','nunique')).reset_index();write_tsv(out,P/'paired_unit_effects.tsv',u);result=[]
 for aset in ef.analysis_set.unique():
  for score in SCORES:
   for A,B in PAIRS:
    name=f'niche{A}_minus_niche{B}'
    for dataset in list(man.dataset_id.unique())+['POOLED_UNITS','DATASET_EQUAL']:
     z=u[u.analysis_set.eq(aset)&u.score.eq(score)&u.pair.eq(name)];z=z[z.dataset_id.eq(dataset)] if dataset not in ['POOLED_UNITS','DATASET_EQUAL'] else z
     # Global collapse true analysis_unit before resampling, even if cross-dataset duplication exists.
     zz=z.groupby('analysis_unit_uid').effect.mean();x=zz.to_numpy();n=len(x);near=1/float(d.AUC_denominator.median());val=z.groupby('dataset_id').effect.mean().to_numpy() if dataset=='DATASET_EQUAL' else x
     ci=np.quantile(rng.choice(x,(2000,n),replace=True).mean(1),[.025,.975]) if n>=5 and dataset!='DATASET_EQUAL' else [np.nan,np.nan]
     result.append(dict(analysis_set=aset,score=score,pair=name,dataset_id=dataset,n_units=n,n_known_units=z.loc[z.known_patient,'analysis_unit_uid'].nunique(),n_anonymous_units=z.loc[~z.known_patient,'analysis_unit_uid'].nunique(),n_datasets=z.dataset_id.nunique(),mean_effect=np.mean(val) if len(val) else np.nan,median_effect=np.median(val) if len(val) else np.nan,min_effect=min(val) if len(val) else np.nan,max_effect=max(val) if len(val) else np.nan,n_positive=int((x>near).sum()),n_negative=int((x< -near).sum()),n_near_zero=int((np.abs(x)<=near).sum()),ci_low=ci[0],ci_high=ci[1],interval_status='UNIT_BOOTSTRAP_EXPLORATORY' if n>=5 and dataset!='DATASET_EQUAL' else 'DESCRIPTIVE_ONLY'))
 rr=pd.DataFrame(result);write_tsv(out,P/'paired_dataset_effects.tsv',rr);write_tsv(out,G/'gene_direct_contrasts.tsv',gene);ll=pd.DataFrame(loo);write_tsv(out,G/'leave_one_gene_out_all20.tsv',ll)
 lu=ll.groupby(['pair','layer','removed_gene','analysis_unit_uid']).agg(original_effect=('original_effect','mean'),renormalized19_effect=('renormalized19_effect','mean'),fixed_denominator_effect=('fixed_denominator_effect','mean')).reset_index();write_tsv(out,G/'leave_one_gene_out_unit_effects.tsv',lu)
 concentration=[]
 for sec,g in d.groupby('section_id'):
  ar=np.load(G/'by_section'/f'{sec}.npz');idx=pd.Index(ar['spot_uid']).get_indexer(g.spot_uid)
  for n,z in g.groupby('accepted_label'):
   ii=pd.Index(ar['spot_uid']).get_indexer(z.spot_uid)
   for layer in ['obs','epi']:
    v=ar[layer+'_contrib'][ii].sum(0);tot=v.sum();concentration.append(dict(section_id=sec,niche_name=f'niche{n}',layer=layer,top_gene=TOP[int(v.argmax())],max_gene_contribution_fraction=v.max()/tot if tot else np.nan,top3_fraction=np.sort(v)[-3:].sum()/tot if tot else np.nan))
 write_tsv(out,G/'contribution_concentration.tsv',concentration)
 # Independent raw rings joined by original UID, never Harmony dimensions.
 orig=pd.read_parquet(OLD/'03_features/common_anchors.parquet',columns=['spot_uid']);f=np.load(OLD/'03_features/features.npz');idx=pd.Index(orig.spot_uid).get_indexer(d.spot_uid);assert (idx>=0).all();cont=[]
 for ring in ['q1','q2']:
  for typ in ['GC B','SPP1 TAM','Neutrophil']:
   d['_ring_value']=f[ring][idx,TYPES.index(typ)]
   for sec,g in d.groupby('section_id'):
    x=g._ring_value.to_numpy();Z=np.column_stack([g.epi_fraction_center_mean,np.log1p(g.model_UMI),g['center_mean_'+typ]]);std=Z.std(0);ok=len(g)>=30 and (x>0).sum()>=20 and x.std()>1e-8 and (std>1e-8).all();D=np.column_stack([np.ones(len(g)),(Z-Z.mean(0))/np.where(std>0,std,1)]);cond=np.linalg.cond(D);ok=ok and cond<1e8
    for score in SCORES[:2]:
     y=g[score].to_numpy();raw=spearmanr(x,y).statistic if x.std()>1e-8 and y.std()>0 else np.nan;adj=np.nan
     if ok:
      xr=rankdata(x);yr=rankdata(y);xr-=D@np.linalg.lstsq(D,xr,rcond=None)[0];yr-=D@np.linalg.lstsq(D,yr,rcond=None)[0];adj=np.corrcoef(xr,yr)[0,1]
     cont.append(dict(section_id=sec,dataset_id=g.dataset_id.iloc[0],analysis_unit_uid=g.analysis_unit_uid.iloc[0],ring=ring,target=typ,score=score,n=len(g),raw_spearman=raw,adjusted_rank_association=adj,design_condition=cond,status='CALCULABLE' if ok else 'LOW_VARIATION_OR_COLLINEAR',interpretation='conditional association; not spatial avoidance/contact'))
 cc=pd.DataFrame(cont);write_tsv(out,P/'continuous_composition_diagnostics.tsv',cc);write_tsv(out,P/'continuous_composition_unit_summary.tsv',cc.groupby(['dataset_id','analysis_unit_uid','ring','target','score'])[['raw_spearman','adjusted_rank_association']].mean().reset_index());write_tsv(out,P/'continuous_composition_dataset_summary.tsv',cc.groupby(['dataset_id','analysis_unit_uid','ring','target','score'])[['raw_spearman','adjusted_rank_association']].mean().reset_index().groupby(['dataset_id','ring','target','score'])[['raw_spearman','adjusted_rank_association']].median().reset_index())
 write_tsv(out,P/'all9_reference_summary.tsv',pd.read_csv(V4/'04_niche_comparison/niche_MP10_summary_all9.tsv',sep='\t'));write_tsv(out,P/'all9_v4_versus_rest_reference.tsv',pd.read_csv(V4/'04_niche_comparison/unit_bootstrap_effects.tsv',sep='\t'))
 for dataset in man.dataset_id.unique():
  write_tsv(out,P/'by_dataset'/dataset/'paired_effects.tsv',rr[rr.dataset_id.eq(dataset)]);fig,axes=plt.subplots(1,3,figsize=(12,4),sharey=True)
  for ax,(A,B) in zip(axes,PAIRS):
   for j,sc in enumerate(SCORES[:2]):
    z=u[u.dataset_id.eq(dataset)&u.analysis_set.eq('MAIN')&u.pair.eq(f'niche{A}_minus_niche{B}')&u.score.eq(sc)];ax.scatter(np.arange(len(z))+.1*j,z.effect,s=20,label='observed' if j==0 else 'model Epi')
   ax.axhline(0,c='gray');ax.set_title(f'niche{A} minus niche{B}');ax.set_xlabel('analysis units');ax.set_ylabel('direct median score difference')
  axes[0].legend();fig.suptitle(dataset);fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(safe(out,out/'06_figures'/dataset/('direct_pairs.'+ext)),transparent=True,dpi=150)
  plt.close(fig)
 safe(out,G/'gene_level_report.md').write_text('Full20 gene counts/ranks/contributions on same per-section model background; formal v4 scores numerically reproduced. B5 excluded from full20 and retained in audit. Attributed nonzero counts are not measured epithelial detection. Leave-one-out includes all20 genes with fixed-background fixed-denominator and renormalized19 diagnostics; none replaces formal scores.\n');log('direct pair comparisons and gene sensitivity complete')
if __name__=='__main__':main(OUT)
