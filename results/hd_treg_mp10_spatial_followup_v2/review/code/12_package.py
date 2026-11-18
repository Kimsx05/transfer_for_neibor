from pathlib import Path
import os,shutil,zipfile,json,subprocess
R=Path(__file__).resolve().parents[1];review=R/'review';review.mkdir(exist_ok=True)
# This directory contains only this task's deliberately selected small deliverables.
for p in review.rglob('*'):
 if p.is_file()or p.is_symlink():p.unlink()
files=list(R.glob('*.md'))+list(R.glob('*.yaml'))+[R/'rerun.sh',R/'replot.sh']+list((R/'code').glob('*.py'))+list((R/'code').glob('*.R'))+list((R/'tables').glob('*.tsv'))+list((R/'figures').glob('*.png'))
parq=['Cervilla_paired_analysis','design_r200_common','design_r200_competition','common_support_reconstructed_centers','common_support_design','expanded_mixtures','mixture_rule_results','original70_replayed']+[s+'_evidence'for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']]
files += [R/'tables'/(s+'.parquet')for s in parq]+[R/'tables/ROI_neighbor_sources.tsv.gz']
for src in files:
 if not src.exists():raise FileNotFoundError(src)
 dst=review/src.relative_to(R);dst.parent.mkdir(parents=True,exist_ok=True)
 try:os.link(src,dst)
 except OSError:shutil.copy2(src,dst)
zipfile_path=R/'review.zip'
def pack():
 with zipfile.ZipFile(zipfile_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6)as z:
  for p in sorted(review.rglob('*')):
   if p.is_file():z.write(p,arcname='hd_treg_mp10_spatial_followup_v2_review/'+str(p.relative_to(review)))
pack();size=int(subprocess.check_output(['du','-sk',str(R)]).split()[0])*1024
text=f'''# 实际交付状态

- 三个任务分支已完成计算与交付；候选校准规则未采用，不等于未执行。
- 24/24固定ROI，512 µm视野，50/200 µm统计半径不变。
- 主共同中心247；计数级共同空间支持170；50 µm共同中心7，仅描述。
- 2700条原混合重放、21600条扩展混合；五片旧/新证据标签完全相同。
- 49项针对性实际验证通过（tables/validation.tsv），Python与R加载器均实际重读。
- 结果目录新增物理磁盘占用约 **{size/1024**2:.1f} MiB**（含review和zip；review文件采用硬链接避免重复占盘）。
- review.zip约 **{zipfile_path.stat().st_size/1024**2:.1f} MiB**。无counts、H5AD、Seurat、表达参考或图像原始副本。
- PNG与小型表可离线审核；加载原表达/重绘H&E依赖manifest中的服务器文件。
- 无正式通讯、细胞身份升级、上游重新分割或新五片连续关联。
- 无排队待运行分支。信息限制：50 µm同中心推断不足，旧完整区域共同支持为空，其他四片未有可复用cell限定MP10，MPBC组织掩膜未知。
'''
(R/'DELIVERY.md').write_text(text);dst=review/'DELIVERY.md'
if dst.exists():dst.unlink()
os.link(R/'DELIVERY.md',dst);pack()
with zipfile.ZipFile(zipfile_path)as z:
 assert z.testzip()is None;assert not any(n.endswith(('.h5ad','.h5','.rds','.RDS','.npz'))for n in z.namelist())
print('PACKAGE_OK',zipfile_path.stat().st_size,'bytes',flush=True)
