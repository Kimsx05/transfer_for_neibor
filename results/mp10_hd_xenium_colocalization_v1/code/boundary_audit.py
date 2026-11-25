from base import *
import gzip,shapely
from display import *
from matplotlib.patches import Polygon
out=[]
for s in XE:
 d=O/'xenium'/s;z=pd.read_parquet(d/'labels.parquet')
 with gzip.open(next((O/'xenium/geo').glob('*_'+s+'_cell_boundaries*')),'rb') as f:b=pd.read_parquet(f)
 b=b.sort_values('cell_id',kind='stable');labs,ids=pd.factorize(b.cell_id,sort=False)
 poly=shapely.polygons(shapely.linearrings(b[['vertex_x','vertex_y']].to_numpy(),indices=labs))
 v=z.loc[ids];pt=shapely.points(v[['x_um','y_um']].to_numpy());dist=shapely.distance(pt,poly);valid=shapely.is_valid(poly);area=shapely.area(poly)
 qc=pd.DataFrame({'polygon_valid':valid,'polygon_area_um2':area,'distance_centroid_to_polygon_um':dist,'polygon_to_original_area':area/v.cell_area.to_numpy(),'cell_area':v.cell_area,'flag':v.flag,'epi_eligible':v.epi_eligible,'chol_group':v.chol_group},index=ids)
 qc['review_boundary']=(~valid)|(dist>2)|(area<=0)
 qc.to_parquet(d/'boundary_geometry.parquet')
 ix=np.array(qc.review_boundary);out.append((s,len(v),(~valid).sum(),(dist>2).sum(),(area<=0).sum(),ix.sum(),np.sum(ix&v.epi_eligible),np.sum(ix&v.chol_group.eq('Chol-high Epi').fillna(False)),float(np.max(dist))))
 if ix.any():
  i=int(np.argmax(np.where(ix,dist,-1)));cx,cy=v[['x_um','y_um']].iloc[i];win=(abs(v.x_um-cx)<=50)&(abs(v.y_um-cy)<=50)
  fig,ax=plt.subplots(figsize=(4,4))
  for k in np.flatnonzero(win):
   ax.add_patch(Polygon(np.asarray(poly[k].exterior.coords),facecolor=colors.get(v.major.iloc[k],'#dddddd'),edgecolor='black' if k==i else 'white',linewidth=1.5 if k==i else .35))
  ax.scatter(cx,cy,c='black',s=8);ax.set(xlim=(cx-50,cx+50),ylim=(cy+50,cy-50),xlabel='x (µm)',ylabel='y (µm)');ax.set_aspect('equal');save(fig,s,'boundary_geometry_review')
put(pd.DataFrame(out,columns=['sample','cells','invalid_polygon','centroid_outside_over_2um','zero_area','boundary_review','eligible_epi_flagged','chol_high_flagged','max_distance_um']),O/'reports/boundary_review.tsv')

