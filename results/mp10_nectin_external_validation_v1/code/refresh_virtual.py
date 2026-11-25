"""Refresh only the optional, read-only interface and the report appendix."""
from pathlib import Path
import csv,json,hashlib,gzip
D=Path(__file__).resolve().parents[1]
V=D.parent/'mp10_srebf_nectin_virtual_perturbation_v1'
F=V/'config/virtual_perturbation_manifest.json'
def run():
 status='pending';rows=[];notes=[];meta={}
 if F.exists():
  meta=json.loads(F.read_text());status=str(meta.get('status',meta.get('completion_status',meta.get('run_status','pending')))).lower()
 if status=='complete':
  for key in ['summary','robustness','coverage']:
   vals=meta.get('files',meta.get('outputs',{}));path=vals.get(key) if isinstance(vals,dict) else None
   entry=next((e for e in vals if isinstance(e,dict) and e.get('path','').endswith('_'+key+'.tsv')),None) if isinstance(vals,list) else None
   if entry:path=entry['path']
   if isinstance(path,dict):path=path.get('path')
   path=path or meta.get(key+'_file') or ('tables/virtual_perturbation_'+key+'.tsv')
   f=(V/path).resolve()
   if not f.is_relative_to(V.resolve()) or not f.is_file():
    notes.append(key+' 缺少有效声明文件');continue
   if entry and hashlib.sha256(f.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Virtual interface fingerprint mismatch '+str(f))
   with (gzip.open(f,'rt') if f.suffix=='.gz' else f.open()) as h:r=list(csv.DictReader(h,delimiter='\t'))
   for x in r:rows.append({'interface_table':key,**x})
   notes.append(key+' '+str(len(r))+' 行')
  if len(notes)!=3 or any('缺少' in s for s in notes):status='incomplete interface'
  else:
   dest=D/'tables/virtual_interface.tsv'
   with dest.open('w') as f:
    w=csv.DictWriter(f,list(dict.fromkeys(k for x in rows for k in x)),delimiter='\t');w.writeheader();w.writerows(rows)
 msg='虚拟接口状态为 '+status+'。'+('；'.join(notes)+'。' if notes else '本次未读取模型输出。')
 if rows:msg+=' 接口按原 assessment_status、metric_type、effect_scale、direction_available、uncertainty_type 保留；网络距离本身无升降方向，模型不确定性与实测干预分别解释。'
 for gene in ['SREBF1','SREBF2']:
  r=next((r for r in rows if r.get('interface_table')=='summary' and r.get('analysis_id')=='main' and r.get('perturbed_gene')==gene and r.get('endpoint_gene')=='NECTIN2'),None)
  if r:
   rr=[v for v in rows if v.get('interface_table')=='robustness' and v.get('perturbed_gene')==gene and v.get('endpoint_gene')=='NECTIN2' and v.get('assessment_status')=='EVALUABLE']
   msg+=f" {gene} 对 NECTIN2 的主模型等级为 {r['assessment_status']}，{r['metric_type']}={float(r['effect_value']):.4g}，单列方向读数为 {r['predicted_direction']}，模型 q={r['q_value']}。{len(rr)} 次可估重拟合有 {sum(v.get('direction_matches_main')=='True' for v in rr)} 次同向；不确定性为 {r['uncertainty_type']}。"
 (D/'qc/virtual_interface.json').write_text(json.dumps(dict(status=status,manifest=str(F),notes=notes),ensure_ascii=False,indent=2))
 p=D/'RESULTS.md'
 if p.exists():
  s=p.read_text();s=s.replace(s[s.index('<!--virtual-->')+14:s.index('<!--/virtual-->')],'\n'+msg+'\n') if '<!--virtual-->' in s and '<!--/virtual-->' in s else s
  p.write_text(s)
 evidence=D/'tables/evidence_by_link.tsv'
 if evidence.exists():
  with evidence.open() as f:items=[r for r in csv.DictReader(f,delimiter='\t') if r['source']!='virtual interface']
  items.append(dict(link='Epi 网络干预预测',source='virtual interface',status=status,result=msg))
  with evidence.open('w') as f:
   w=csv.DictWriter(f,['link','source','status','result'],delimiter='\t');w.writeheader();w.writerows(items)
 return msg
if __name__=='__main__':print(run())

