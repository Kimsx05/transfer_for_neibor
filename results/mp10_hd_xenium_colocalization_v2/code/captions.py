from common import *
import matplotlib.pyplot as plt
from figures import save
def captions():
 lines=['# 图件阅读说明','每张图为单独 PNG/PDF，Provisional / Codex-reviewed，人工待审。横向 x 向右、纵向 y 向下，空间单位均为 µm。']
 records=[]
 for s in HD+XE:
  hd=s in HD;v,d=folders(s);out=OUT/'figures'/s;rows=pd.read_csv(d/'score_QC.tsv',sep='\t')
  unit='Treg component bins/mm²' if hd else 'Treg cells/mm²'
  labs={'density':unit,'Treg/T':'Treg / reliable T','Treg/all':'Treg / all neighbors','any contact':'fraction with boundary proximity','restricted distance um':'restricted distance (µm)'}
  eff=pd.read_csv(d/'effects.tsv.gz',sep='\t')
  for end,lab in labs.items():
   sub=eff[(eff.endpoint==end)&(eff.block_um==400)&(eff.identity_definition=='supported')]
   if sub.empty:continue
   tick=sub.state_definition+' · '+sub.unit_resolution.astype(str)+' · '+sub.radius_or_contact_um.astype(str)+' µm';fig,ax=plt.subplots(figsize=(5.7,max(2.5,len(sub)*.38+1.3)));yy=np.arange(len(sub));ax.hlines(yy,sub.CI_low,sub.CI_high,color='#547384');ax.scatter(sub.effect,yy,s=18,color='#842d50');ax.axvline(0,ls='--',color='#999',lw=.7);ax.set(yticks=yy,yticklabels=tick,xlabel='High − low, '+lab);save(fig,out/('effect_'+end.replace('/','_').replace(' ','_')))
  lines.append('## '+s+'\n\n[逐图入口](../figures/'+s+'/index.html)。')
  for r in rows.itertuples():
   method='历史 AUCell 口径，全 18,085 feature，K905' if hd else ('冻结 H05 参考支持度' if r.state=='H05_transfer' else 'UCell，共同 5,001 基因，maxRank 1,500')
   lines.append(f'- {r.res} · {r.state} 使用{method}。本片 Q95={r.q95:.6g}，严格大于阈值的 {r.high:,} 个为 high，其余 {r.low:,} 个合格上皮为 low。')
   records.append((s,r.res,r.state,method,r.q95,r.high,r.low,r.QC))
  lines.append('邻域主读数按 50 µm、逐中心等权计算。区间为 400 µm 空间块重采样的 95% 区间。Xenium 近接触报告可用边界；HD 密度单位为支持 bin/mm²。散点的平滑线是等频分箱均值的平滑趋势，niche 编号在片内解释。')
 (OUT/'reports/FIGURES.md').write_text('\n\n'.join(lines)+'\n');put(pd.DataFrame(records,columns=['sample','resolution','program','method','Q95','high','low','QC']),OUT/'reports/figure_thresholds.tsv')
if __name__=='__main__':captions()
