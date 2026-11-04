from report import *
def make(sec):
 rr=M.set_index('section_id').loc[sec];cached=O/'figures/sections'/rr.dataset_id/(sec+'_4_full_range_scatter.pdf')
 if all(cached.with_name(cached.stem+suf).exists() for suf in ['.pdf','.png','_panel_only.pdf','_legend_only.pdf']):return dict(section_id=sec,dataset_id=rr.dataset_id,page='supplement_full_range',pdf=str(cached),png=str(cached.with_suffix('.png')))
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');a=pd.read_parquet(O/'tables'/f'{sec}_pearson.parquet');de=pd.read_parquet(O/'tables'/f'{sec}_delta.parquet');r=M.set_index('section_id').loc[sec];normal=r.source_context=='adjacent_normal_source';sel=pd.read_csv(O/'reports'/f'{sec}_scatter_selection.tsv',sep='\t');f,axs=plt.subplots(2,2,figsize=(18,13));notes=[]
 for i,rr in enumerate(sel.itertuples()):notes.append(scatter_pair(axs[i],d,a,de,rr.candidate_type,rr.program,normal,range_type='full_range'))
 f.suptitle(sec+' | raw spot Pearson | full observed candidate-distance range (no clipping)',fontsize=14);f.text(.06,.025,'\n'.join(notes)+'\nCommon-support formal comparison is separate in the 3-page main report. Pearson P does not correct spatial dependence.',fontsize=11);leg=f.add_axes([.4,.002,.3,.02]);leg.axis('off');leg.legend(handles=[Patch(color=COL[0],label='P95'),Patch(color=COL[1],label='nonP95'),Patch(color='#ff9900',label='distance 0')],ncol=3,frameon=False);f.subplots_adjust(top=.91,bottom=.17,hspace=.36,wspace=.25);fp=export(f,'sections/'+r.dataset_id+'/'+sec+'_4_full_range_scatter',[leg]);return dict(section_id=sec,dataset_id=r.dataset_id,page='supplement_full_range',pdf=fp,png=fp[:-4]+'.png')
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as ex:save(list(ex.map(make,M.section_id)),O/'tables/full_range_scatter_index.tsv')
