# Reuse imports, fixed palette and transparent PDF/PNG helpers without rerunning all plots.
exec(open('mp10_subtype_screen_round1/scripts/05_abundance_figures.py').read().split('comps=')[0])
sel=pd.read_csv(T/'selected_axes.tsv',sep='\t');su=pd.read_csv(T/'candidate_union_sample_expression_support.tsv.gz',sep='\t');r=pd.read_csv(T/'proportion_models.tsv',sep='\t')
dscols={d:to_hex(plt.get_cmap('tab10')(i)) for i,d in enumerate(sorted(s.dataset.unique()))}
for typ in sel.candidate_subtype.unique():
 q=p[p.level.eq('subtype')&p.subtype.eq(typ)].copy();q=q.merge(s[['sample_uid','MP10_eligible']],on='sample_uid');q['MP10_eligible']=q.MP10_eligible.astype(str).str.lower().eq('true')
 for ds in ['ALL']+sorted(s.dataset.unique()):
  z=q if ds=='ALL' else q[q.dataset.eq(ds)]
  groups=sorted(z.dataset.unique())
  for comp,(col,values) in {'Normal_Tumor':('tissue_status',['Normal','Tumor']),'T1_T2':('T_category',['T1','T2']),'MP10_low_high':('MP10_group',['low','high'])}.items():
   a=z[z[col].isin(values)]
   if comp=='T1_T2':a=a[a.tissue_status.eq('Tumor')]
   if a.empty:continue
   fig,axs=plt.subplots(len(groups),2,figsize=(10,max(3.7,3.1*len(groups))),squeeze=False)
   for i,d in enumerate(groups):
    for j,den in enumerate(['parent','TME']):
     b=a[a.dataset.eq(d)&a.denominator_type.eq(den)&a.proportion.notna()];ax=axs[i,j]
     for pat,rr in b.groupby('patient_uid'):
      if rr[col].nunique()==2:
       rr=rr.sort_values(col,key=lambda x:x.map(dict(zip(values,[0,1]))));ax.plot(rr[col].map(dict(zip(values,[0,1]))),rr.proportion*100,c='#999999',alpha=.5,lw=.7,zorder=1)
     for k,g in enumerate(values):
      bb=b[b[col].eq(g)].sort_values('sample_uid');jit=np.linspace(-.12,.12,len(bb));ax.scatter(k+jit,bb.proportion*100,s=30,c=dscols[d],edgecolors='white',linewidths=.3,zorder=3)
     ax.set(xticks=[0,1],xticklabels=values,xlim=(-.45,1.45),ylabel=f'{typ} / {den} (%)',title=d);ax.grid(axis='y',alpha=.15)
   fig.suptitle(f'{typ}: {comp}\nOne point per specimen; lines use recorded patient mappings',y=1.01);fig.tight_layout();save(fig,ds,'candidate_points_'+safe(typ)+'_'+comp)
  a=z[z.MP10_eligible&z.proportion.notna()]
  if not a.empty:
   fig,axs=plt.subplots(len(groups),2,figsize=(10,max(3.7,3.1*len(groups))),squeeze=False)
   for i,d in enumerate(groups):
    for j,den in enumerate(['parent','TME']):
     b=a[a.dataset.eq(d)&a.denominator_type.eq(den)];ax=axs[i,j];ax.scatter(b.MP10_sample,b.proportion*100,c=dscols[d],s=32,edgecolors='white',linewidths=.3)
     if len(b)>=3 and b.MP10_sample.nunique()>1:
      xx=np.linspace(b.MP10_sample.min(),b.MP10_sample.max(),30);fit=np.polyfit(b.MP10_sample,b.proportion*100,1);ax.plot(xx,np.polyval(fit,xx),c=dscols[d],lw=1)
     ax.set(xlabel='Original sample MP10 mean',ylabel=f'{typ} / {den} (%)',title=f'{d}; n={len(b)}');ax.grid(alpha=.15)
   fig.suptitle(f'{typ}: continuous MP10\nWithin-dataset descriptive trends; adjusted inference in model table',y=1.01);fig.tight_layout();save(fig,ds,'candidate_continuous_'+safe(typ))
for _,a in sel.iterrows():
 q=su[su.LR_ID.eq(a.LR_ID)&su.source.eq(a.source)&su.target.eq(a.target)&su.MP10_group.isin(['low','high'])]
 for ds in ['ALL']+sorted(s.dataset.unique()):
  z=q if ds=='ALL' else q[q.dataset.eq(ds)];groups=sorted(z.dataset.unique())
  if not groups:continue
  fig,axs=plt.subplots(len(groups),2,figsize=(10,max(3.5,2.8*len(groups))),squeeze=False)
  for i,d in enumerate(groups):
   b=z[z.dataset.eq(d)]
   for j,(val,det,role) in enumerate([('ligand_mean','ligand_min_detection',f'{a.ligand} in sender'),('receptor_mean','receptor_min_detection',f'{a.receptor} in receiver')]):
    ax=axs[i,j]
    for k,g in enumerate(['low','high']):
     bb=b[b.MP10_group.eq(g)&b.evaluable].sort_values('sample_uid');xx=k+np.linspace(-.13,.13,len(bb));ax.scatter(xx,bb[val],s=20+bb[det]*100,c=dscols[d],alpha=.8,edgecolors='black',linewidths=.3)
    ax.set(xticks=[0,1],xticklabels=['low','high'],xlim=(-.45,1.45),ylabel='Mean log-normalized RNA',title=f'{d}: {role}');ax.grid(axis='y',alpha=.15)
    if not b.evaluable.any():ax.text(.5,.5,'Source/target insufficient',transform=ax.transAxes,ha='center')
  fig.suptitle(f'{a.LR_display}: {a.source} → {a.target}\nEach point: evaluable specimen; size: minimum subunit detection fraction',y=1.015);fig.tight_layout();save(fig,ds,'sample_support_'+safe(a.candidate_subtype)+'_'+a.LR_ID)
print('Candidate plots complete')
