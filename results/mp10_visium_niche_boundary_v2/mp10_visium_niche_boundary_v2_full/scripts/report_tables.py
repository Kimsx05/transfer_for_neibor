from lib import *
import argparse
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');v=p.parse_args().version;out=O if v=='primary' else O/'versions'/v
f=out/'tables/group_effects.tsv.gz';df=read(f);df['q_two_boundaries']=np.nan
for _,idx in df[df.scope.eq('overall')].groupby(['variant','width','endpoint','group']).groups.items():df.loc[idx,'q_two_boundaries']=bh(df.loc[idx,'p'],2)
save(df,f)
cl=read(out/'classification/section_classification.tsv')
rows=[]
for r in cl.itertuples():
    rows.append((r.mode,r.provisional_architecture_group,r.section_id,r.group_reason,r.confidence,r.scale1000))
doc=['# 逐切片临时标签','', '主规则为 500 μm。最后一列列出 1000 μm 敏感性标签。正常切片单独投影。','', '|队列|主标签|切片及审图|判定理由|置信度|1000 μm|','|---|---|---|---|---|---|']
for md,g,s,why,cf,k in sorted(rows):doc.append(f'|{md}|{g}|[{s}](../review/{s}.jpg)|{why}|{cf}|{k}|')
(out/'classification/SAMPLES.md').write_text('\n'.join(doc)+'\n')

