from common import *
from report import merged,read,put,md,bh
def run():
 d=OUT/'reports';pri=read(d/'primary.tsv');a=pri[pri.endpoint=='density'][['sample','state_definition','effect','FDR']].rename(columns={'state_definition':'state','effect':'density_delta','FDR':'density_q'})
 sq=merged('squidpy_state.tsv')
 if len(sq):
  a=a.merge(sq[sq.group=='high'][['sample','state','Z','OE']],how='left')
 c=merged('CRAWDAD_formal.tsv')
 if len(c):
  cc=c.pivot(index=['sample','state'],columns='tile_um',values='mean_proportion_Z');cc.columns=['CRAWDAD_Z_'+str(r) for r in cc];a=a.merge(cc.reset_index(),how='left')
 fixed=merged('fixed_target_curve.tsv')
 if len(fixed):
  f=fixed[fixed.r==50].pivot(index=['sample','state'],columns='group',values='ratio_to_unit_prevalence').reset_index().rename(columns={'high':'high_edge_ratio_to_prevalence','low':'low_edge_ratio_to_prevalence'});a=a.merge(f,how='left')
 curves=read(d/'distance_global.tsv')
 if len(curves):
  cv=curves[(curves.definition=='supported')&(curves.block==400)&(curves.res.astype(str)!='32')].drop_duplicates(['sample','state']);a=a.merge(cv[['sample','state','global_P','FDR']].rename(columns={'FDR':'curve_q'}),how='left')
 native=[];influence=[];gt=[];bounds=[]
 for s in HD+XE:
  v,ds=folders(s);L=read(ds/'leave_one.tsv.gz')
  if len(L):
   for kk,t in L.groupby(['sample','res','state','endpoint']):influence.append((*kk,t.full.iloc[0],t.leave_one.min(),t.leave_one.max(),t.loc[(t.leave_one-t.full).abs().idxmax(),'excluded_id'],((t.leave_one*t.full)<0).any()))
  for key in (['MP10'] if s in HD else STATES):
   en=read(ds/f'spatstat_envelope_{key}.tsv')
   if len(en):native.append((s,key,len(en),((en.obs<en.lo)|(en.obs>en.hi)).any()))
  if s in XE:
   q=pd.read_parquet(v/'boundary_geometry.parquet');p=pd.read_parquet(ds/'boundary_pairs.parquet');bounds.append((s,len(q),q.review_boundary.sum(),(~q.polygon_valid).sum(),p.gap_um.eq(0).sum(),(p.gap_um<=1).sum(),(p.gap_um<=2).sum(),len(p)))
  e=read(ds/'effects.tsv.gz');src=read(ds/'score_QC.tsv')
  for _,row in src.iterrows():
   met=e[(e.unit_resolution.astype(str)==str(row.res))&(e.state_definition==row.state)&(e.endpoint=='density')&(e.identity_definition=='supported')&(e.radius_or_contact_um==50)&(e.block_um==400)]
   if len(met):gt.append((s,row.res,row.state,row.n,row.high,met.effect.iloc[0],met.balanced_effect.iloc[0],met.Spearman.iloc[0],row.QC))
 if native:a=a.merge(pd.DataFrame(native,columns=['sample','state','curve_radii','spatstat_outside_global_band']),how='left')
 put(a,d/'tool_comparison.tsv');put(pd.DataFrame(influence,columns=['sample','res','state','endpoint','full','LOO_min','LOO_max','most_influential_target','sign_flip']),d/'leave_one_summary.tsv')
 put(pd.DataFrame(bounds,columns=['sample','polygons','review_flag','invalid','zero_gap_edges','edges_1','edges_2','edges_5']),d/'boundary_audit.tsv')
 put(pd.DataFrame(gt,columns=['sample','resolution','state','eligible','high','density_effect','balanced_effect','Spearman','score_QC']),d/'sensitivity.tsv')
if __name__=='__main__':run()
