from base import *
from display import *
from scipy.stats import mannwhitneyu
from matplotlib.patches import Polygon
import gzip
def review(s):
 d=O/'xenium'/s;z=pd.read_parquet(d/'labels.parquet');x,ids,g=tenx(next((P/'highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5').glob('*_'+s+'_*h5')));assert np.array_equal(ids,z.index)
 idx=g.feature_type.eq('Gene Expression').to_numpy();x=x[:,idx];genes=pd.Index(g.name[idx]);maps=json.loads((O/'shared/major_markers.json').read_text());markers=set(sum(maps.values(),[]))
 tx=pd.read_csv(O/'shared/taxonomy.tsv',sep='\t');markers.update(' '.join(tx.markers).split());markers.update('KRT20 UPK1A UPK2 UPK3A GATA3 FOXA1 TP63 KRT5 KRT14 KRT17 DSG3 IVL KRT6A'.split())
 markers=sorted(markers&set(genes));vals=pd.DataFrame(x[:,genes.get_indexer(markers)].toarray(),index=ids,columns=markers)
 clusters=[]
 for cl,b in z.groupby('cluster'):
  a=vals.loc[b.index];freq=b.major_candidate.value_counts(normalize=True)
  for gene in a:
   clusters.append((cl,len(b),freq.index[0],freq.iloc[0],b.flag.eq('subtype_supported').mean(),gene,(a[gene]>0).mean(),a[gene].mean()))
 put(pd.DataFrame(clusters,columns=['cluster','n','leading_major_candidate','candidate_fraction','subtype_supported_fraction','gene','detection','mean_raw']),d/'cluster_marker_review.tsv.gz')
 e=z[z.epi_eligible].copy();r=e[['Chol_biosynthesis_8','H05_transfer']].rank(pct=True);e['discordance']=np.select([(r.H05_transfer>=.8)&(r.Chol_biosynthesis_8<=.2),(r.H05_transfer<=.2)&(r.Chol_biosynthesis_8>=.8)],['high_transfer_low_measured','low_transfer_high_measured'],default='remainder')
 cols=[g for g in 'KRT20 UPK1A UPK2 UPK3A GATA3 FOXA1 TP63 KRT5 KRT14 KRT17 DSG3 IVL KRT6A'.split() if g in vals]
 expr=e[['x_um','y_um','nGeneTranscripts','nFeature','cell_area','H05_transfer','Chol_biosynthesis_8','MP10_covered','major_support','discordance','cluster']].join(vals[cols]);expr[expr.discordance!='remainder'].to_parquet(d/'discordant_epithelial_identity.parquet')
 rows=[]
 for k,a in expr.groupby('discordance'):
  for gene in cols:rows.append((k,len(a),gene,(a[gene]>0).mean(),a[gene].mean(),a.nGeneTranscripts.median(),a.cell_area.median(),a.major_support.median()))
 put(pd.DataFrame(rows,columns=['group','n','epi_marker','marker_fraction','mean_raw','median_transcripts','median_area','median_major_support']),d/'discordance_review.tsv')
 # Representative polygons selected by area and mixture evidence, independent of cholesterol and MP10.
 with gzip.open(next((O/'xenium/geo').glob('*_'+s+'_cell_boundaries*')),'rb') as f:b=pd.read_parquet(f)
 xy=z[['x_um','y_um']].to_numpy()
 candidates={'large_area':z.cell_area.idxmax(),'small_area':z.cell_area.idxmin()}
 mix=z[z.flag=='mixed']
 if len(mix):candidates['mixed']=mix.sort_values('major_runner',ascending=False).index[0]
 roi=[]
 for role,cid in candidates.items():
  center=z.loc[cid,['x_um','y_um']].to_numpy(float);j=np.flatnonzero((np.abs(xy-center)<=50).all(1));v=z.iloc[j];poly=b[b.cell_id.isin(v.index)]
  fig,ax=plt.subplots(figsize=(4,4)); 
  for c,a in poly.groupby('cell_id'):
   ax.add_patch(Polygon(a[['vertex_x','vertex_y']].to_numpy(),facecolor=colors.get(v.loc[c,'major'],'#aaa'),edgecolor='black' if c==cid else 'white',linewidth=1.5 if c==cid else .35))
  ax.set_xlim(center[0]-50,center[0]+50);ax.set_ylim(center[1]+50,center[1]-50);ax.set_aspect('equal');ax.set(xlabel='x (µm)',ylabel='y (µm)');save(fig,s,'segmentation_'+role)
  roi.append((role,cid,*center,z.loc[cid,'cell_area'],z.loc[cid,'flag'],len(j)))
 put(pd.DataFrame(roi,columns=['selection','center_id','x_um','y_um','area_um2','flag','cells_in_window']),d/'segmentation_ROIs.tsv')
 # Shared 80-um grid retains raw data at discordant locations and locally comparable continuous scores.
 e['col']=np.floor(e.x_um/80);e['row']=np.floor(e.y_um/80);grid=e.groupby(['row','col']).agg(n=('H05_transfer','size'),H05_transfer=('H05_transfer','mean'),Chol_biosynthesis_8=('Chol_biosynthesis_8','mean'),MP10_covered=('MP10_covered','mean'))
 genes2=json.loads((O/'shared/programs.json').read_text());raw=pd.DataFrame(x[:,genes.get_indexer(sorted(set(sum(genes2.values(),[]))))].toarray(),index=ids,columns=sorted(set(sum(genes2.values(),[]))))
 grid=grid.join(raw.loc[e.index].groupby([e.row,e.col]).sum().add_prefix('raw_'))
 ok=grid.n>=5;grid['rank_gap']=np.nan;grid.loc[ok,'rank_gap']=grid.loc[ok,'H05_transfer'].rank(pct=True)-grid.loc[ok,'Chol_biosynthesis_8'].rank(pct=True)
 grid.to_parquet(d/'local_state_expression_80um.parquet')
if __name__=='__main__':
 for s in sys.argv[1:]:review(s)

