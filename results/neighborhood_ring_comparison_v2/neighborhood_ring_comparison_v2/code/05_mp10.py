from common import *
from itertools import combinations

def main():
 f=OUT/'05_comparison/frozen_manifest.json';review=OUT/'01_regions/LATEST_REVIEW.json'
 if not f.exists() or not review.exists():log('M06 BLOCKED: blind labels and actual region review must be frozen; no spatial MP10 values read');return
 freeze=json.loads(f.read_text());rpath=Path(json.loads(review.read_text())['path']);r=json.loads(rpath.read_text())
 if not freeze['blind_labels_frozen'] or not r['MP10_unblinding_allowed']:log('M06 BLOCKED: region review incomplete; no spatial MP10 values read');return
 # Validated versioned review only; original feature/cache hashes remain unchanged.
 m=pd.read_parquet(rpath.parent/'reviewed_master.parquet');m=m[m.common_anchor].reset_index(drop=True)
 mf=pd.read_csv(ROOT/'neighborhood_v1/00_inputs_and_audit/metadata/frozen_MP_scores_manifest.csv');sc=[]
 for row in mf.itertuples():
  t=pd.read_csv(row.score_file,usecols=['barcode','EpiMP10_Top20_AUCell']);t['section_id']=row.section;sc.append(t)
 m=m.merge(pd.concat(sc),on=['section_id','barcode'],how='left',validate='1:1');assert m.EpiMP10_Top20_AUCell.notna().all();av=pd.read_csv(OUT/'00_preflight/MP10_gene_availability.tsv',sep='\t');av['gene_stratum']=av.n_available.astype(str)+'::'+av.missing_genes.fillna('');m=m.merge(av[['section','gene_stratum']],left_on='section_id',right_on='section',validate='m:1');high=json.loads((OUT/'00_preflight/input_contract.json').read_text())['high_epi_threshold'];all_effects=[];descs=[]
 for name in ['R1','R2','R12']:
  for alg in ['KM','LEIDEN']:
   config=name+'_'+alg;l=pd.read_parquet(OUT/'04_clustering'/name/alg/'labels.parquet');t=m.merge(l[['spot_uid','cluster_label_native_application','cluster_label_common_transfer']],on='spot_uid',validate='1:1')
   strata={'all':np.ones(len(t),bool),'tumor_supported':t.region_final.isin(['tumor_supported_epithelium','tumor_supported_region']),'normal_epithelium':t.region_final.eq('normal_epithelium_supported'),'reviewed_unknown':t.region_final.eq('uncertain_or_conflicting'),'complete_rings':t.complete_rings,'high_epi':t.epithelial_fraction>=high,'reliable_depth_same_subset':t.depth_valid}
   for lname in ['cluster_label_native_application','cluster_label_common_transfer']:
    for st,mask in strata.items():
     tt=t[mask].copy()
     for (ds,s,c),g in tt.groupby(['dataset_id','section_id',lname]):descs.append(dict(config=config,label=lname,stratum=st,dataset_id=ds,section_id=s,cluster=c,n=len(g),mean_MP10=g.EpiMP10_Top20_AUCell.mean(),median_MP10=g.EpiMP10_Top20_AUCell.median(),patient_uid=g.patient_uid.iloc[0]))
     for (ds,s),g in tt.groupby(['dataset_id','section_id']):
      for ca,cb in combinations(sorted(t[lname].unique()),2):
       h=g[g[lname].isin([ca,cb])].dropna(subset=['epithelial_fraction','nCount']);sizes=h[lname].value_counts()
       if len(sizes)!=2 or sizes.min()<30:continue
       z=(h[lname].to_numpy()==cb).astype(float);y=h.EpiMP10_Top20_AUCell.to_numpy();epi=h.epithelial_fraction.to_numpy();umi=np.log1p(h.nCount.to_numpy());designs={'unadjusted':np.c_[np.ones(len(h)),z],'epi_UMI_adjusted':np.c_[np.ones(len(h)),z,epi,umi]}
       if st=='reliable_depth_same_subset':designs['epi_UMI_depth_adjusted']=np.c_[np.ones(len(h)),z,epi,umi,h.depth_hops.to_numpy()]
       for model,X in designs.items():
        if np.linalg.matrix_rank(X)<X.shape[1]:continue
        b=np.linalg.lstsq(X,y,rcond=None)[0];all_effects.append(dict(config=config,label=lname,stratum=st,model=model,dataset_id=ds,section_id=s,patient_uid=h.patient_uid.iloc[0],source_patient_identity_known=bool(h.source_patient_identity_known.iloc[0]),gene_stratum=h.gene_stratum.iloc[0],cluster_a=ca,cluster_b=cb,n=len(h),effect_b_minus_a=b[1]))
 tab(OUT/'06_mp10/section_descriptions.tsv',descs);ef=pd.DataFrame(all_effects);tab(OUT/'06_mp10/section_effects.tsv',ef)
 if len(ef):
  known=ef[ef.source_patient_identity_known & ef.patient_uid.notna()];keys=['config','label','stratum','model','dataset_id','gene_stratum','cluster_a','cluster_b'];pat=known.groupby(keys+['patient_uid'],dropna=False).effect_b_minus_a.mean().reset_index();tab(OUT/'06_mp10/patient_effects.tsv',pat);summ=[];rng=np.random.default_rng(123)
  for key,g in pat.groupby(keys,dropna=False):
   vals=g.effect_b_minus_a.to_numpy();row=dict(zip(keys,key));row.update(n_patients=len(vals),mean_effect=vals.mean(),CI_low=np.nan,CI_high=np.nan,status='descriptive_less_than_5_patients')
   if len(vals)>=5:
    b=np.array([rng.choice(vals,len(vals),replace=True).mean() for _ in range(2000)]);row.update(CI_low=np.quantile(b,.025),CI_high=np.quantile(b,.975),status='patient_bootstrap_pointwise_CI_exploratory_no_Pvalues')
   summ.append(row)
  tab(OUT/'06_mp10/dataset_patient_effects.tsv',summ)
 js(OUT/'06_mp10/STATUS.json',dict(status='EXPLORATORY_PATIENT_EFFECTS_COMPLETE',region_version=r['version'],inference='Within-section point estimates; equal sections within true patient; patient bootstrap within dataset+gene availability; no spot-level standard errors/Pvalues; no cross-dataset pooled score',multiple_comparisons='pointwise exploratory intervals, no confirmatory significance claims; no P values',causal_claims=False,remaining='source validation and matched-depth composition-only controls remain mandatory'))
if __name__=='__main__':main()
