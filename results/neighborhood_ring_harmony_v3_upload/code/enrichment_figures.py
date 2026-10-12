from base import *
import matplotlib.pyplot as plt
F=OUT/'04_comparison/figures'
for p in sorted((OUT/'03_main').iterdir()):
 pr=pd.read_csv(p/'profiles_R1_R2.tsv',sep='\t');fig,axes=plt.subplots(2,1,figsize=(14,max(7,pr.cluster.nunique()*.4)))
 for ax,ring in zip(axes,['R1','R2']):
  a=pr[pr.ring==ring].pivot(index='cluster',columns='cell_type',values='log2_enrichment').reindex(columns=TYPES);im=ax.imshow(a,aspect='auto',cmap='RdBu_r',vmin=-3,vmax=3);ax.set_xticks(range(30),TYPES,rotation=70,ha='right',fontsize=6);ax.set_yticks(range(len(a)),a.index.astype(str));ax.set_title(p.name+' '+ring+' relative enrichment (original composition)');fig.colorbar(im,ax=ax,shrink=.6,label='log2((niche+1e-6)/(global+1e-6))')
 fig.tight_layout();fig.savefig(F/(p.name+'__enrichment.png'),dpi=150,transparent=True);fig.savefig(F/(p.name+'__enrichment.pdf'),transparent=True);plt.close(fig)
log('18 enrichment figure pairs complete')
