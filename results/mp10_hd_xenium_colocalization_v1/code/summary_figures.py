from base import *
from display import *
def plot_results(s):
 hd=s in HD;d=O/('hd' if hd else 'xenium')/s;z=pd.read_parquet(d/'labels.parquet')
 if not hd:
  cs=z.chol_group.copy();cs=cs.fillna('Other tissue');cs[z.major=='Unknown']='Unknown';cs[z.major=='Mixed']='Mixed';colors['Other tissue']='#dedede'
  spatial(s,z,'Chol_high_Epi',cs)
 for r in [50,200]:
  t=pd.read_csv(d/f'effects_{r}.tsv',sep='\t').set_index('neighbor').reindex(['Treg','Endothelial (all)','Fibroblast (all)','Myeloid (all)'])
  fig,ax=plt.subplots(figsize=(5.1,3.3));pos=np.arange(len(t));h=t.fraction_high*100;o=t.fraction_other*100
  for i in pos:
   if np.isfinite(o.iloc[i]) and np.isfinite(h.iloc[i]):ax.plot([o.iloc[i],h.iloc[i]],[i,i],color='#777777',lw=1.1)
  ax.scatter(o,pos,s=25,color='#ce9968',label='other Epi');ax.scatter(h,pos,s=30,color='#a61d48',label=('H05' if hd else 'Chol-high Epi'))
  ax.set_yticks(pos,['Treg','Endothelial','Fibroblast','Myeloid']);ax.invert_yaxis();ax.set_xlabel(f'Mean neighbor fraction at {r} µm (%)');ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1));save(fig,s,f'neighborhood_composition_{r}')
 fig,ax=plt.subplots(figsize=(5,3.4))
 for k in ['subtype_supported','major_only','uncertain','mixed']:
  v=z.loc[z.flag==k,'major_support'].dropna().sort_values()
  if len(v):ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',label=k.replace('_',' '),color={'subtype_supported':'#228833','major_only':'#4477aa','uncertain':'#888888','mixed':'#ee7733'}[k])
 ax.legend(frameon=False,fontsize=8);ax.set_xlabel('Leading major weight' if hd else 'Leading major transfer support');ax.set_ylabel('Cumulative fraction');ax.set_xlim(0,1);ax.set_ylim(0,1);save(fig,s,'support_distribution')
if __name__=='__main__':
 for s in sys.argv[1:]:plot_results(s)

