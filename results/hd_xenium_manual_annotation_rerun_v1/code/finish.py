"""Verify all feature blocks, primary IDs and manual templates before publishing the workbook."""
from data import *
from importlib import metadata
from openpyxl import Workbook, load_workbook
from openpyxl.styles import PatternFill,Font,Alignment
from openpyxl.utils import get_column_letter
from scipy.stats import false_discovery_control
import platform

manual='manual_major manual_subtype manual_label supporting_markers contradictory_markers confidence reviewer notes'.split()
audit=[];inventory=[];manifests=[];priorities={};w=Workbook();w.remove(w.active)
for p,ss in S.items():
 d=O/p;mf=pd.read_csv(d/'input_manifest.tsv',sep='\t');manifests.append(mf)
 x,g,src=combine(p);n,ng=x.shape
 obj=pd.read_parquet(d/'objects.parquet').set_index('cell_key',verify_integrity=True)
 assert set(obj.index)==set(src.cell_key)
 obj=obj.loc[src.cell_key]
 for field in ['sample','original_cell_id','source_row','n_counts','nFeature','x_um','y_um','source']:
  assert np.array_equal(obj[field].to_numpy(),src[field].to_numpy())
 for s in ss:
  allq=pd.read_parquet(d/(s+'_qc.parquet'));keep=allq[allq.qc_pass]
  assert allq.cell_key.is_unique and np.array_equal(allq.source_row,np.arange(len(allq)))
  assert np.array_equal(keep.cell_key,obj.index[obj['sample']==s])
  assert keep.shape[0]==mf.loc[mf['sample']==s,'included'].item()
 audit.append((p,'raw matrix, full QC, composite IDs and source rows',n,'passed'))
 lat=pd.read_parquet(d/'latent.parquet').loc[obj.index]
 um=obj[['UMAP1','UMAP2']].to_numpy()
 assert np.isfinite(um).all() and np.linalg.matrix_rank(um-um.mean(0))==2 and not np.allclose(um,lat.iloc[:,:2])
 graph=sparse.load_npz(d/'expression_graph.npz')
 rows=pd.read_csv(d/'graph_rows.tsv.gz',sep='\t');assert rows.cell_key.tolist()==lat.index.tolist()
 assert graph.shape==(n,n) and not graph.diagonal().any() and (graph!=graph.T).nnz==0
 cfg=json.loads((d/'graph_config.json').read_text())
 assert cfg['n_edges']==graph.nnz//2 and cfg['seed']==1729 and cfg['embedding'].endswith('simplicial_set_embedding')
 del graph,lat
 audit.append((p,'UMAP, weighted symmetric expression graph, node index',cfg['n_edges'],'passed'))
 for r in R:
  dr=d/('r'+r);col=cname(r);cc=obj[col]
  labels=pd.read_parquet(dr/'labels.parquet').loc[obj.index];assert np.array_equal(cc,labels[col])
  tally=pd.crosstab(cc,obj['sample']).reindex(columns=ss,fill_value=0)
  comp=pd.read_csv(dr/'cluster_composition.tsv',sep='\t').set_index('cluster')
  assert np.array_equal(tally.to_numpy(),comp[[s+'_n' for s in ss]].to_numpy()) and comp.n_objects.sum()==n
  colors=pd.read_csv(dr/'colors.tsv',sep='\t');assert np.array_equal(colors.cluster,comp.index) and colors.color.is_unique
  sizes=comp.n_objects; pooled_n=within_n=ntest=nz=bad=small=0
  feat=g.feature_id.to_numpy()
  for kind in ['pooled','within_sample']:
   cols=['cluster','feature_id','n_in','n_out','mean_norm_in','mean_norm_out','pct_in','pct_out','log2FC','status']
   if kind=='pooled':cols+=['p_value','p_adj_bh']
   else:cols+=['sample']
   got=set()
   for frame in pd.read_csv(dr/(kind+'_all_genes.tsv.gz'),sep='\t',usecols=cols,chunksize=ng):
    assert np.array_equal(frame.feature_id,feat) and frame.cluster.nunique()==1
    c=int(frame.cluster.iat[0]);a=frame.n_in.iat[0];b=frame.n_out.iat[0]
    if kind=='pooled':
     key=c;assert a==sizes.loc[c] and a+b==n
     ok=frame.status.eq('ok')
     assert frame.loc[~ok,['p_value','p_adj_bh']].isna().all().all()
     assert frame.loc[ok,'p_value'].between(0,1).all() and frame.loc[ok,'p_adj_bh'].between(0,1).all()
     if ok.any():assert np.allclose(false_discovery_control(frame.loc[ok,'p_value']),frame.loc[ok,'p_adj_bh'],rtol=2e-9,atol=1e-300)
     ntest+=ok.sum();nz+=frame.status.eq('all_zero').sum();bad+=(~ok).sum();pooled_n+=len(frame)
    else:
     s=frame['sample'].iat[0];key=(c,s)
     assert frame['sample'].nunique()==1 and a==tally.loc[c,s] and b==tally[s].sum()-a
     assert (frame.status.str.startswith('insufficient_')==(min(a,b)<20)).all()
     if min(a,b)<20:
      assert frame[['mean_norm_in','mean_norm_out','pct_in','pct_out','log2FC']].isna().all().all();small+=1
     within_n+=len(frame)
    assert key not in got;got.add(key)
    if kind=='pooled' or min(a,b)>=20:
     assert frame.pct_in.between(0,1).all() and frame.pct_out.between(0,1).all()
     assert np.allclose(np.log2((frame.mean_norm_in+1e-9)/(frame.mean_norm_out+1e-9)),frame.log2FC,atol=1e-9)
    if kind=='pooled' and c==0:
     aa=(cc==c).to_numpy();jj=np.linspace(0,ng-1,9,dtype=int);meas=x[:,jj].toarray()
     real=meas*10000./src.n_counts.to_numpy()[:,None]
     assert np.allclose(real[aa].mean(0),frame.mean_norm_in.iloc[jj])
     assert np.allclose((meas[~aa]>0).mean(0),frame.pct_out.iloc[jj])
   assert len(got)==len(sizes)*(1 if kind=='pooled' else len(ss))
  inventory.append(dict(platform=p,resolution=r,n_clusters=len(sizes),n_genes=ng,pooled_rows=pooled_n,within_sample_rows=within_n,tested_rows=ntest,all_zero_rows=nz,untestable_rows=bad,within_insufficient_cluster_samples=small))
  audit.append((p,'full feature IDs, all cluster/sample blocks, effects, BH r'+r,pooled_n+within_n,'passed'))
 template=pd.read_csv(d/'manual_review.tsv',sep='\t')
 assert template[manual].isna().all().all() and len(template)==obj.leiden_r0_5.nunique()
 sheet=w.create_sheet(p.upper() if p=='hd' else 'Xenium');sheet.append(template.columns.tolist())
 for row in template.itertuples(index=False,name=None):sheet.append([None if pd.isna(v) else v for v in row])
 for i,cell in enumerate(sheet[1],1):
  cell.fill=PatternFill('solid',fgColor='2D674F' if cell.value in manual else '284C69');cell.font=Font(color='FFFFFF',bold=True)
  sheet.column_dimensions[get_column_letter(i)].width=65 if str(cell.value).startswith('Top') else max(14,len(str(cell.value))+3)
 sheet.freeze_panes='D2';sheet.auto_filter.ref=sheet.dimensions
 for row in sheet.iter_rows(min_row=2):
  for cell in row:
   if str(sheet.cell(1,cell.column).value).startswith('Top'):cell.alignment=Alignment(wrap_text=True,vertical='top')
 sheet.sheet_format.defaultRowHeight=44
 priorities[p]=pd.read_csv(d/'review_priorities.tsv',sep='\t')
 ff=O/'figures'/p;assert {y.stem for y in ff.glob('*.png')}=={y.stem for y in ff.glob('*.pdf')}
 for nm in ['UMAP_sample']+['UMAP_r'+r for r in R]+[s+'_spatial_r0.5' for s in ss]+[s+'_UMAP_r0.5' for s in ss]:assert (ff/(nm+'.pdf')).stat().st_size>1000
 for pat in ['positive_dotplot','positive_heatmap','Treg_markers']:assert list(ff.glob(pat+'*.pdf'))
 audit.append((p,'manual columns blank and figure pairs present',len(list(ff.glob('*.png'))),'passed'))
 del x,obj
w.save(O/'manual_annotation.xlsx')
for sh in load_workbook(O/'manual_annotation.xlsx',read_only=True).worksheets:
 rowit=sh.values;header=next(rowit);assert all(row[header.index(k)] is None for row in rowit for k in manual)
mf=pd.concat(manifests,ignore_index=True);write(mf,O/'input_manifest.tsv')
mi=pd.DataFrame(inventory);write(mi,O/'tables/deg_inventory.tsv');write(pd.DataFrame(audit,columns=['platform','check','n','status']),O/'tables/validation.tsv')
qcols=['platform','sample','objects','included','excluded','low_features','low_counts','base_geometry_flag','secondary_image_suspect','area_outlier','boundary_degenerate'];write(mf[qcols],O/'tables/qc_summary.tsv')
versions={pkg:metadata.version(pkg) for pkg in ['scvi-tools','torch','lightning','anndata','numpy','scipy','pandas','pyarrow','igraph','leidenalg','umap-learn','numba','matplotlib','openpyxl']}
versions.update(python=platform.python_version(),python_executable=sys.executable);(O/'tables/versions.json').write_text(json.dumps(versions,indent=2))
def mdtable(f):
 return '\n'.join(['| '+' | '.join(map(str,f.columns))+' |','| '+' | '.join(['---']*len(f.columns))+' |']+['| '+' | '.join(map(str,row))+' |' for row in f.itertuples(index=False,name=None)])
report=['# 表达聚类与人工注释准备','', '主版本为 Leiden 0.5。打开 [人工注释工作簿](manual_annotation.xlsx)，按平台与 cluster 对照 [图件](figures/) 和逐样本效应。绿色表头下的人工字段留空。','',mdtable(mf[qcols]),'', '过滤原因可重叠。逐对象原始行号、过滤原因及重建来源存于各平台的样本 QC 文件。','',mdtable(mi),'','## 主分辨率的复核线索','']
for p,t in priorities.items():
 report+=['### '+p,'']
 for key,name in [('dominant_sample_fraction','单一样本占比'),('two_epi_and_two_immune','同一对象内多项上皮与免疫 marker 共同检出'),('median_nFeature','表达复杂度')]:
  sort=t.sort_values(key,ascending=key=='median_nFeature').head(3)
  report += [name+'：'+'，'.join(f"簇 {int(v.cluster)}（{getattr(v,key):.3g}，n={int(v.n_objects)}）" for v in sort.itertuples())+'。','']
 template=pd.read_csv(O/p/'manual_review.tsv',sep='\t')
 ranked=t.nsmallest(2,'median_nFeature').merge(template[['cluster','Top10_positive']],on='cluster',validate='one_to_one')
 for v in ranked.itertuples():
  report += [f"簇 {int(v.cluster)} 的中位检出 feature 数为 {v.median_nFeature:.1f}，正 marker 摘要为 {v.Top10_positive}。",'']
 report+=['结合各样本在平台总体中的基线占比解读，稀有簇须检查样本内两侧是否均达到 20 个对象。共同检出可由混合、组织背景或真实表达造成，仅作为定位线索。','']
report+=['## 模型与图的实际规模','']
for p in S:
 d=O/p
 tr=pd.read_csv(d/'train_loss.tsv',sep='\t');va=pd.read_csv(d/'validation_loss.tsv',sep='\t')
 graph=json.loads((d/'graph_config.json').read_text())
 report+=[f"{p} 训练 {len(tr)} 个 epoch，末轮 train loss {tr.iloc[-1,-1]:.3f}、validation loss {va.iloc[-1,-1]:.3f}。表达图含 {graph['n_edges']:,} 条无向边。完整轨迹存于平台目录。",'']
report+=['逐簇详情见各平台 review_priorities.tsv、r0.5/cluster_composition.tsv 和样本内全基因表。统计用于聚类后人工筛选 marker；同一切片内对象存在依赖，不能据此推出患者层面的重复性。','']
(O/'REPORT.md').write_text('\n'.join(report))
print('verified both platforms and workbook',flush=True)

