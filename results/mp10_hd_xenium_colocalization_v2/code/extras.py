from common import *
from scipy.spatial.distance import jensenshannon
def run():
 review=[];quality=[];mix=[];states=[];fit=[];interp=[];panrows=[]
 for s in HD+XE:
  v,d=folders(s);hd=s in HD;z=labels(s);idt=pd.read_parquet(d/'identity.parquet');p=OUT/'niche'/s;n=pd.read_parquet(p/'MENDER.parquet').niche
  e=pd.read_parquet(d/('centers_16.parquet' if hd else 'centers_cell.parquet'))
  genes=(v/'counts/genes.txt').read_text().splitlines() if hd else (V/'shared/common_panel.txt').read_text().splitlines()
  cc=pd.read_parquet(paths(s)/'candidates.parquet')
  if hd:
   ww=pd.read_parquet(v/'weights.parquet').reindex(z.index);g=sparse.load_npz(v/'graph_50.npz');fitted=ww.notna().all(axis=1).to_numpy('i4');den=g@fitted;mat=g@ww.fillna(0).to_numpy('f4');mat/=np.where(den>0,den,np.nan)[:,None]
   agg=pd.DataFrame(mat,columns=ww.columns,index=z.index).groupby(n).mean()
   for k,row in agg.iterrows():
    labs=sorted(row.index,key=lambda x:-row[x]);fit.append((s,k,n.eq(k).sum(),np.mean(den[n.eq(k)]>0),labs[0],row[labs[0]],row['Treg'],row.Epithelial))
    for lab,x in row.items():mix.append((s,k,lab,x))
  else:
   for key in STATES:
    for h in [True,False]:
     sel=e[key+'_high'].eq(h);ids=e.index[sel]
     for meas in STATES[:2]:
      n0=e.loc[ids,meas+'_detected'];fr=e.loc[ids,meas+'_dominance']
      states.append((s,key,'high' if h else 'low',meas,len(ids),n0.eq(0).mean(),n0.eq(1).mean(),n0.ge(2).mean(),fr.mean()))
  candidates=idt[idt.candidate_union]
  for mode,zz in candidates.groupby(['grade','review_mode']):review.append((s,*mode,len(zz)))
  for k,nn in n.groupby(n):
   ix=nn.index;sel=e.index.isin(ix)
   quality.append((s,k,len(ix),idt.UMI.loc[ix].median(),idt.loc[ix].major.isin(['Unknown','Mixed']).mean(),len(e[sel]),e.loc[sel,'coverage_50'].mean(),idt.grade.loc[ix].eq('supported').sum()))
  for name in ['MENDER','BuildNicheAssay','CellCharter_all','CellCharter_masked']:
   kfile=p/(name+('.tsv.gz' if name=='BuildNicheAssay' else '.parquet'));status=kfile.exists()
   if not status:continue
   m=pd.read_csv(kfile,sep='\t').set_index('id') if name=='BuildNicheAssay' else pd.read_parquet(kfile)
   assert m.index.is_unique and set(m.index)==set(z.index)
  for k in CORE+PAN:panrows.append((s,k,k in genes))
 put(pd.DataFrame(fit,columns=['sample','niche','nodes','fitted_neighborhood_fraction','leading_weight','mean_leading','mean_Treg','mean_Epi']),OUT/'reports/niche_HD_weight_summary.tsv')
 put(pd.DataFrame(mix,columns=['sample','niche','subtype','mean_neighbor_RCTD_weight']),OUT/'reports/niche_HD_weights.tsv')
 put(pd.DataFrame(states,columns=['sample','classification','group','measured_program','cells','zero_members','one_member','multigene','mean_max_share']),OUT/'reports/high_program_support.tsv')
 put(pd.DataFrame(quality,columns=['sample','niche','nodes','median_UMI','unknown_mixed_fraction','epi_centers','mean_coverage','Treg_supported']),OUT/'reports/niche_quality.tsv')
 put(pd.DataFrame(review,columns=['sample','grade','review_mode','n']),OUT/'reports/review_modes.tsv')
if __name__=='__main__':run()
