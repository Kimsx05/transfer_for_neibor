from run import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from PIL import Image
import subprocess,gzip,json
font_manager.fontManager.addfont(subprocess.check_output(['fc-match','-f','%{file}','Arial'],text=True))
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'axes.linewidth':1.1,'lines.linewidth':1.1,'figure.facecolor':'white'})
def save(fig,name,axes,legends=[]):
 p=O/'figures'/name
 fig.savefig(str(p)+'.pdf',facecolor='white',bbox_inches='tight');fig.savefig(str(p)+'.png',dpi=300,facecolor='white',bbox_inches='tight')
 if legends:
  for a in legends:a.set_visible(False)
  fig.savefig(str(p)+'_panel_only.pdf',bbox_inches='tight',facecolor='white')
  for a in legends:a.set_visible(True)
  for a in axes:a.set_visible(False)
  for t in fig.texts:t.set_visible(False)
  fig.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight',facecolor='white')
 plt.close(fig)
def heat(mat,title,name,cap=.3):
 fig,ax=plt.subplots(figsize=(max(8,len(mat.columns)*1.2),max(5,len(mat)*.27+2)))
 m=ax.pcolormesh(np.ma.masked_invalid(mat.to_numpy(float)),vmin=-cap,vmax=cap,cmap='RdBu_r')
 ax.set_xticks(np.arange(len(mat.columns))+.5,mat.columns,rotation=60,ha='right');ax.set_yticks(np.arange(len(mat))+.5,[s.replace('HALLMARK_','').replace('REACTOME_','') for s in mat.index]);ax.invert_yaxis();ax.set_title(title);ax.set_facecolor('#ddd')
 cb=fig.colorbar(m,ax=ax);cb.solids.set_rasterized(False);fig.tight_layout();save(fig,name,[ax],[cb.ax])
def summary():
 a=pd.read_csv(O/'tables/summary.tsv',sep='\t');a=a[(a.source.str.startswith('tumor'))&(a.population=='units')]
 for kind,pr,name in [('distance',PROG,'program_distance'),('interaction',PROG,'interactions')]:
  z=a[(a.kind==kind)&a.program.isin(pr)&a.axis.str.endswith('P95')].pivot(index='program',columns='axis',values='effect')
  heat(z,'P95 / adjusted effects / equal datasets of analysis units',name,.3 if kind=='distance' else .08)
 z=a[(a.kind=='distance')&a.program.eq(RAW)].set_index('axis')[['effect','low','high']];heat(z,'MP10 adjusted distance correlations and hierarchical 95% intervals','MP10_effects',.2)
 z=a[a.axis.eq('MP10')].set_index('program')[['effect','low','high']];heat(z,'MP10 versus all 55 spot programs','MP10_programs',.5)
 z=a[(a.kind=='distance')&a.axis.isin(['local_interface','observed_outer','strict_outer'])].pivot(index='program',columns='axis',values='effect');heat(z,'Local interface and lesion-edge references','interface_programs',.3)
 t=pd.read_csv(O/'tables/targets.tsv.gz',sep='\t');z=t[t['quantile']==95].pivot(index='section_id',columns='target',values='abundance_median');heat(np.log1p(z),'P95 target estimated abundance (log1p median)','target_support',2)
def atlas(sec):
 r=M.set_index('section_id').loc[sec];d=pd.read_parquet(O/'spots'/f'{sec}.parquet')
 fig,aa=plt.subplots(3,3,figsize=(15,16));aa=aa.ravel()
 for ax,t in zip(aa,list(TARGETS)):
  ax.scatter(d.x,d.y,s=3,c='#ddd');mask=d[t+'_P95_target'];ax.scatter(d.loc[mask,'x'],d.loc[mask,'y'],s=5,c='#b33e65');ax.set_title(t+' P95');ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off')
 fig.suptitle(sec);save(fig,sec+'_targets',list(aa))
 fig,aa=plt.subplots(1,3,figsize=(17,6));aa[0].imshow(Image.open(r.image_path));aa[0].set_title('H&E / cyan local, orange lesion edge')
 with gzip.open(ROOT/'results/mp10_pathology_regions_v3_4_final/by_dataset'/r.dataset_id/sec/'geometry.json.gz','rt') as h: geo=json.load(h)
 outer=pd.DataFrame(geo['roi_edges']['data'],columns=geo['roi_edges']['columns']);outer=outer[outer.boundary_type=='candidate_observed_outer']
 if len(outer):aa[0].add_collection(LineCollection(outer[['x1_px','y1_px','x2_px','y2_px']].to_numpy().reshape(-1,2,2),colors='#D55E00',linewidths=1.2))
 edges=pd.read_csv(O/'tables'/f'{sec}_edges.tsv',sep='\t') if (O/'tables'/f'{sec}_edges.tsv').stat().st_size>2 else pd.DataFrame()
 if len(edges):
  lines=[]
  # μm to image mapping estimated from existing corresponding coordinates.
  affine=np.linalg.lstsq(np.c_[d.x_um,d.y_um,np.ones(len(d))],d[['x','y']],rcond=None)[0]
  for row in edges.to_numpy():
   lines.append(np.c_[np.array(row[2:6],float).reshape(2,2),np.ones(2)]@affine)
  aa[0].add_collection(LineCollection(lines,colors='#08b5b5',linewidths=1.2))
 for ax,col,title in [(aa[1],'local_side','Supported sides (+1 epi, -1 stroma)'),(aa[2],'local_interface_um','Signed local distance (μm)')]:
  ax.scatter(d.x,d.y,s=3,c='#ddd');ok=d[col].notna();im=ax.scatter(d.loc[ok,'x'],d.loc[ok,'y'],s=4,c=d.loc[ok,col],cmap='coolwarm');ax.set_title(title);ax.invert_yaxis()
  cb=fig.colorbar(im,ax=ax);cb.solids.set_rasterized(False)
 for ax in aa:ax.set_aspect('equal');ax.axis('off')
 fig.suptitle(sec);save(fig,sec+'_interface',list(aa),[a for a in fig.axes if a not in aa])
 return sec
if __name__=='__main__':
 summary()
 with ProcessPoolExecutor(max_workers=6) as p:
  for s in p.map(atlas,M.section_id):print(s,flush=True)
