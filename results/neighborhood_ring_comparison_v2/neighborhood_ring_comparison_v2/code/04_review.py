from common import *
import argparse
from collections import deque
ALLOWED=['tumor_supported_epithelium','normal_epithelium_supported','non_epithelial_region','uncertain_or_conflicting','tumor_supported_region','stroma_supported','other_supported']

def distances(g,allowed,seeds):
 d=np.full(g.shape[0],np.nan);queue=deque(int(i) for i in np.flatnonzero(seeds));d[seeds]=0
 while queue:
  i=queue.popleft()
  for j in g.indices[g.indptr[i]:g.indptr[i+1]]:
   if allowed[j] and np.isnan(d[j]):d[j]=d[i]+1;queue.append(int(j))
 return d

def apply_review(path):
 base=pd.read_csv(OUT/'01_regions/ROI_review.tsv',sep='\t');r=pd.read_csv(path,sep='\t').fillna('');assert r.roi_id.is_unique
 assert set(r.roi_id).issubset(set(base.roi_id)), 'Unknown ROI IDs'
 reviewed=r[r.reviewed_label.ne('')].copy();assert reviewed.reviewed_label.isin(ALLOWED).all(),'Invalid label'
 assert reviewed[['reviewer','review_date','review_comment']].ne('').all().all(),'Actual reviewer, date and evidence comment required'
 pd.to_datetime(reviewed.review_date,errors='raise');version='review_'+sha(path)[:12];od=OUT/'01_regions'/version;od.mkdir(exist_ok=True)
 idx=reviewed.set_index('roi_id');allm=[]
 for ds in sorted(base.dataset_id.unique()):
  for p in sorted((OUT/'01_regions'/ds).glob('*__regions.parquet')):
   m=pd.read_parquet(p);sel=m.roi_id.isin(idx.index);m.loc[sel,'region_final']=m.loc[sel,'roi_id'].map(idx.reviewed_label);m.loc[sel,'region_review_status']='reviewed';m['region_version']=version
   s=m.section_id.iloc[0];g=sparse.load_npz(OUT/'02_graphs'/ds/(s+'__R1.npz'));tum=m.region_final.isin(['tumor_supported_epithelium','tumor_supported_region']).to_numpy();non=m.region_final.isin(['normal_epithelium_supported','non_epithelial_region','stroma_supported','other_supported']).to_numpy();known=tum|non
   inner=tum&((g@non.astype(int))>0);outer=non&((g@tum.astype(int))>0)
   di=distances(g,tum,inner);do=distances(g,non,outer);depth=np.where(tum,di,do)
   # Censor distances when an unknown region or physical truncation is no farther than a known interface.
   trunc=(np.asarray(g.sum(1)).ravel()<6)|((g@(~known).astype(int))>0)
   dt=distances(g,known,trunc&known);valid=known&np.isfinite(depth)&(~np.isfinite(dt)|(depth<dt))
   m['inner_boundary']=inner;m['outer_boundary']=outer;m['depth_hops']=np.where(valid,depth,np.nan);m['depth_valid']=valid;m['boundary_status']=np.where(valid,'reviewed_known_interface','unknown_or_truncated');m['boundary_type']='unknown'
   for i in np.flatnonzero(inner|outer):
    ni=g.indices[g.indptr[i]:g.indptr[i+1]];labs=set(m.region_final.iloc[ni]) if inner[i] else {m.region_final.iloc[i]};types=[]
    if 'stroma_supported' in labs:types.append('tumor-stroma')
    if 'normal_epithelium_supported' in labs:types.append('tumor-normal_epithelium')
    if labs.intersection({'other_supported','non_epithelial_region'}):types.append('tumor-other')
    if types:m.loc[i,'boundary_type']=';'.join(types)
   m.to_parquet(od/(s+'.parquet'),index=False);allm.append(m)
 master=pd.concat(allm,ignore_index=True);common=master[master.common_anchor].copy();all_reviewed=master.region_review_status.eq('reviewed').all();master.to_parquet(od/'reviewed_master.parquet',index=False)
 js(od/'review_manifest.json',dict(version=version,source=str(Path(path).resolve()),source_sha256=sha(path),reviewed_ROIs=len(reviewed),pending_ROIs=len(base)-len(reviewed),all_regions_reviewed=bool(all_reviewed),MP10_unblinding_allowed=bool(all_reviewed),anchors_unchanged=True,labels_refit=False,depth_definition='known inner and outer boundary each side=0; hops within side; censored if unknown/gap no farther'))
 js(OUT/'01_regions/LATEST_REVIEW.json',dict(path=str(od/'review_manifest.json')));log(f'review applied: {len(reviewed)}/{len(base)}; labels unchanged')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--review',required=True);a=p.parse_args();apply_review(a.review)
