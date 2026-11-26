from common import *
from spatial import groups
for s in sys.argv[1:] or HD+XE:
 v,d=folders(s);z=labels(s);i=pd.read_parquet(d/'identity.parquet');res=16 if s in HD else 'cell';e=pd.read_parquet(d/f'centers_{res}.parquet');m=['MP10'] if s in HD else STATES
 n=pd.read_parquet(d/f'neighborhood_{res}_200.parquet').set_index('id').loc[e.index];n50=pd.read_parquet(d/f'neighborhood_{res}_50.parquet').set_index('id').loc[e.index];r=e[['x_um','y_um']].copy();r['stratum']=groups(e,400);r['cf200']=n.valid;r['distance']=n.supported_distance_um.fillna(201);r['density']=np.where(n50.valid,n50.supported/e.area_50*1e6,np.nan);r['Tfrac']=np.where(n50.valid,n50.supported/n50.reliable_T.replace(0,np.nan),np.nan)
 t=z[['x_um','y_um']].copy();t['Treg']=i.grade.eq('supported');t['major']=z.major;t['known']=~z.major.isin(['Mixed','Unknown'])
 for key in m:r[key]=e[key+'_high'].astype(int);t[key]=r[key].reindex(t.index).fillna(-1).astype(int)
 put(r.reset_index(),d/'spatstat_input.tsv.gz');put(t.reset_index(),d/'crawdad_input.tsv.gz')
