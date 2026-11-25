from base import *
def read(p):return pd.read_csv(p,sep='\t')
def interpret():
    d=O/'reports';z=read(d/'xenium_states.tsv');v=read(d/'heldout_validation.tsv')
    identity=read(d/'identity_coverage.tsv');eff=read(d/'per_sample_state_effects.tsv')
    notes=['# 首轮判读','', '## 细类支持和测量分母','']
    focus=identity[identity.flag=='subtype_supported']
    for plat,ss in [('Xenium',XE),('HD',HD)]:
        done=[s for s in ss if s in set(identity['sample'])]
        by=focus[focus['sample'].isin(done)].pivot_table(index='supported',columns='sample',values='units',aggfunc='sum',fill_value=0).reindex(columns=done,fill_value=0)
        robust=by.index[(by>=10).sum(axis=1)==len(ss)].tolist() if len(done)==len(ss) else []
        notes.append(f'{plat} 已取得 {len(done)}/{len(ss)} 个样本的完整身份表。'+('每片至少 10 个模型及 marker 双重支持单元的亚群为 '+ '、'.join(robust)+'。' if robust else '覆盖分母见逐片表。'))
        if 'Treg' in by.index:notes.append('Treg 支持单元数为 '+ '，'.join(f'{s} {int(by.loc["Treg",s])}' for s in done)+'。')
    rows=[]
    for s in XE:
        p=O/'xenium'/s/'labels.parquet'
        if not p.exists():continue
        a=pd.read_parquet(p);c=a[a.epi_eligible];high=c.chol_group=='Chol-high Epi';flags=(a.flag=='uncertain')|(a.major=='Unknown')
        q=read(O/'xenium'/s/'discordance_review.tsv') if (O/'xenium'/s/'discordance_review.tsv').exists() else pd.DataFrame()
        for group,f in [('all supported Epi',np.ones(len(c),bool)),('Chol-high Epi',high)]:
            b=c[f];rows.append((s,group,len(b),flags.mean(),(b.Chol_biosynthesis_8_detected==0).mean(),(b.Chol_biosynthesis_8_detected==1).mean(),(b.Chol_biosynthesis_8_detected>=2).mean(),b.Chol_biosynthesis_8_max_gene_fraction.median(),b.nGeneTranscripts.median(),b.cell_area.median()))
    ev=pd.DataFrame(rows,columns=['sample','group','n','all_cell_uncertain_fraction','zero_program_fraction','single_gene_fraction','multi_gene_fraction','median_dominant_gene_share_nonzero','median_transcripts','median_area'])
    put(ev,d/'real_epithelial_score_support.tsv')
    notes+=['','## 三种上皮读数','',f'实际五片内，合成–MP10_covered 的 Spearman ρ 为 {z["合成–MP10 ρ"].min():.3f} 至 {z["合成–MP10 ρ"].max():.3f}，合成–H05_transfer 为 {z["合成–H05 ρ"].min():.3f} 至 {z["合成–H05 ρ"].max():.3f}。前两个集合共享 HMGCR、MVD、FDPS、FDFT1、SQLE 五个基因。合成读数与总 transcripts 也同向（ρ {z["合成–深度 ρ"].min():.3f} 至 {z["合成–深度 ρ"].max():.3f}），局部组织密度和深度需要和分数一起解释。']
    a=v[v['mode']=='thinned'].groupby('score').average_precision.median();notes.append(f'整来源留出加真实深度降采样后，区分原 H05 的 AP 中位数依次为 MP10_covered {a.MP10_covered:.3f}、8 基因合成 {a.Chol_biosynthesis_8:.3f} 和 H05_transfer {a.H05_transfer:.3f}，原始阳性率中位数 {v.prevalence.median():.3f}。保留实测程序作为空间主读数。H05_transfer 适合作为有不确定度的连续对照。')
    notes+=['','## 状态周围的样本差异','']
    for plat,ss,r in [('Xenium',XE,50),('HD',HD,200)]:
        ep=eff[(eff.platform==plat)&(eff.radius==r)]
        if ep.empty:continue
        for k in ['Treg','Endothelial (all)','Fibroblast (all)','Myeloid (all)']:
            t=ep[(ep.neighbor==k)&(ep.section_nodes>=10)&(ep.n_high>=5)&(ep.n_other>=5)].dropna(subset=['delta'])
            if len(t):
                desc='，'.join(f'{b.sample} {100*b.delta:+.3f}' for b in t.itertuples());notes.append(f'{plat} {r} µm 的 {k} 邻域份额差（百分点）为 {desc}，{len(t)}/{len(ss)} 片有足够描述覆盖。')
    pair=read(d/'cross_platform_subtype_pairs.tsv')
    uniform=pair[(pair.cross=='方向一致')&(pair.hd_n==3)&(pair.xenium_n==5)]
    if len(uniform):
        examples=uniform[uniform.hd=='富集'].head(4)
        notes+=['','50 µm 全部八片均可评价且方向相同的亚群对共 '+str(len(uniform))+' 对。'+('其中同向富集有 '+'、'.join(t.a+' 与 '+t.b for t in examples.itertuples())+'。' if len(examples) else '')+'这些是相对于全片标签置换的方向，不评价物理接触或患者间效应。']
    elif set(HD)<=set(identity['sample']):notes+=['','未得到在全部三片 HD 和五片 Xenium 都达到支持门槛且同号的 50 µm 亚群对。较少切片可评价的关系单列于方向表。']
    rel=pd.read_csv(d/'sample_subtype_pairs.tsv.gz',sep='\t')
    pair_ex=rel[(rel['from']=='C1QC TAM')&(rel.to=='Matrix_CAF')&(rel.radius==50)]
    if len(pair_ex)==8:
        vals=pair_ex[pair_ex.platform=='Xenium'].observed_expected
        h=pair_ex[pair_ex.platform=='HD'].set_index('sample').reindex(HD)
        notes.append(f'50 µm 的 C1QC TAM–Matrix_CAF 在五份 Xenium 的观察/置换期望比范围为 {vals.min():.2f}–{vals.max():.2f}；三张 HD 按 Cervilla、D1_Johnson、UC_SCC 顺序分别为 '+ '、'.join(f'{x:.2f}' for x in h.observed_expected)+'。这是 Xenium 内同向的描述性关系，跨平台有异质性。')
    he=eff[(eff.platform=='HD')&(eff.radius==200)&(eff.neighbor=='Matrix_CAF')].set_index('sample')
    if len(he)==3:
        notes.append('200 µm 内冻结 H05 相对其他合格 Epi64 的 Matrix_CAF 邻居份额差为 '+ '，'.join(f'{k} {100*he.loc[k,"delta"]:+.2f} 个百分点' for k in HD)+'。三片均已剔除共享表达来源配对。')
    notes+=['','## 下一步优先级','','首先补 Treg 与相近活化 T 细胞、低覆盖 EC/DC 的身份验证，并核对候选区原分割与原 marker。其次在同区域补充胆固醇程序的直接测量，重点选 H05_transfer 与实测程序秩差大的区域；优先参考三状态空间图和 80 µm 格的原计数。只有在身份分母、状态、密度及跨切片方向均支持时，再选定具体细胞对做独立组织验证。']
    up=d/'unknown_marker_review.tsv'
    if up.exists():
        u=read(up);u=u[u.candidate_major=='Epithelial']
        notes.append('补注释可优先复核强上皮候选但缺少主 marker 的单元。在当前 Unknown 中，候选支持 ≥0.8、优势 ≥0.3，且另有至少两枚既有上皮状态 marker 的细胞数为 '+ '，'.join(f'{x["sample"]} {int(x["two_extra_markers"])}' for _,x in u.iterrows())+'。这些候选均保留在审阅表。')
    sp=d/'reference_source_review.tsv'
    if sp.exists():
        u=read(sp);w=u.pivot_table(index=['sample','reference_source'],columns='discordance',values='mean_source_weight')
        if 'high_transfer_low_measured' in w and 'remainder' in w:
            diff=(w.high_transfer_low_measured-w.remainder).dropna()
            if len(diff):
                (sample,source)=diff.abs().idxmax()
                notes.append(f'对高 H05 支持而低实测合成的细胞，复用相同上皮参考与 panel 投影转移参考来源标签。与其余细胞相比，最大来源平均权重差为 {sample} 对 {source} 的 {diff.loc[(sample,source)]:+.3f}。这个诊断提示参考来源结构的可能影响；其大小与正负方向逐样本报告，不能据此断定配对或来源导致状态。')
    (d/'INTERPRETATION.md').write_text('\n\n'.join(x for x in notes if x)+'\n')
if __name__=='__main__':interpret()

