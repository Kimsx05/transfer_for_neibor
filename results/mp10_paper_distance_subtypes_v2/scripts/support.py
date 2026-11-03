from plots import *
def main():
 rows=[];nr=[]
 for sec in M.section_id:
  try:
   z=pd.read_csv(O/'tables'/f'{sec}_curves.tsv.gz',sep='\t')
   d=pd.read_parquet(O/'spots'/f'{sec}.parquet');mid=[]
   for axis,g in z.groupby('axis'):
    col={'local_interface':'local_interface_um','observed_outer':'distance_to_observed_outer_candidate_um','strict_outer':'signed_distance_um'}.get(axis,axis+'_distance')
    b=pd.cut(d[col],[-np.inf,-1000,-500,-200,0,100,200,500,1000,np.inf]).astype(str)
    for label in g.bin:mid.append(dict(axis=axis,bin=label,distance_median=d.loc[d['query']&b.eq(label),col].median()))
   z=z.merge(pd.DataFrame(mid),on=['axis','bin']);rows.append(z)
   z=pd.read_csv(O/'tables'/f'{sec}_niche_interaction.tsv.gz',sep='\t');nr.append(z)
  except pd.errors.EmptyDataError:pass
 c=pd.concat(rows).merge(M[['section_id','dataset_id','source_context']],on='section_id')
 c=c[c.source_context.str.startswith('tumor')]
 # Equal section mean then equal dataset mean. No pooled spot association.
 g=c.groupby(['dataset_id','axis','bin'])[PROG+[RAW,'distance_median']].mean().groupby(['axis','bin']).mean()
 axes=[t+'_P95' for t in list(TARGETS)[:7]]+['local_interface']
 for page in range(7):
  modules=PROG[page*8:(page+1)*8];fig,axs=plt.subplots(4,2,figsize=(15,18));axs=axs.ravel()
  for ax,pr in zip(axs,modules):
   for axis in axes:
    if axis not in g.index.get_level_values(0):continue
    z=g.loc[axis].copy()
    def order(b):return float(b.split(',')[0].strip('[('))
    z=z.loc[sorted(z.index,key=order)]
    ax.plot(z.distance_median,z[pr],marker='.',label=axis)
   ax.set_title(pr.replace('HALLMARK_','').replace('REACTOME_',''),fontsize=11);ax.set_ylabel('Epi-query median UCell')
   ax.set_xlabel('Mean within-bin query distance (μm)')
  for ax in axs[len(modules):]:ax.set_visible(False)
  leg=fig.add_axes([.12,.01,.78,.05]);leg.axis('off');handles,labels=axs[0].get_legend_handles_labels();leg.legend(handles,labels,ncol=4,loc='center',fontsize=10)
  fig.tight_layout(rect=[0,.08,1,1]);save(fig,f'curves_{page+1}',list(axs),[leg])
 if nr:
  n=pd.concat(nr).merge(M[['section_id','dataset_id']],on='section_id');n['scoped_niche']=np.where(n.niche_version.eq('local_niche_id'),n.dataset_id+' / '+n.niche.astype(str),n.niche.astype(str));write(n,O/'tables/niche_interactions.tsv.gz')
  z=n.groupby(['niche_version','scoped_niche','axis']).beta.mean().unstack('axis')
  heat(z.set_axis([str(a)+' / '+str(b) for a,b in z.index]),'Within-niche mean interaction coefficient (descriptive)','niche_interactions',.1)
  # Keep dataset-local labels scoped by section/dataset in raw tables.
 n=pd.read_csv(O/'tables/niche.tsv.gz',sep='\t')
 z=n[(n.axis=='local_interface')&(n.version=='R12')].merge(M[['section_id','dataset_id','source_context']],on='section_id')
 z=z[z.source_context.str.startswith('tumor')]
 z=z.groupby(['dataset_id','bin','niche']).fraction.mean().groupby(['bin','niche']).mean().unstack()
 if len(z):heat(z,'R12 composition along signed local interface','interface_niche',1)
 genes=set(Path(O/'sc_genes.txt').read_text().splitlines()) if (O/'sc_genes.txt').exists() else set()
 write([dict(program=p,n_original=g.gene_symbol.nunique(),n_matched=len(set(g.gene_symbol)&genes)) for p,g in SETS.groupby('gs_name')],O/'tables/sc_gene_coverage.tsv')
if __name__=='__main__':main()
