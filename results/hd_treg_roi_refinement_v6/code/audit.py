from base6 import *
rows=[];fields=[]
for s in sections:
 d=pd.read_parquet(V5/s/'cells_scored.parquet');g=d.in_D & ~d.image_review_exclusion & ~d.geometry_review & ~d.image_confirmed_invalid & (d.nuclear_bins>0)
 p=g & ~d.mixed & d.pred.eq('Treg') & d.pTreg.ge(.5) & d.core_width.ge(1)
 f=g & d.T_width.ge(2) & d.FOXP3.gt(0)
 for k,v in {'wide2':p&d.T_width.ge(2),'wide1':p&d.T_width.ge(1),'wide1_CD2_only':p&d.T_width.eq(1)&d.CD2.gt(0),'FOXP3_before_mixed':f,'FOXP3_after_mixed':f&~d.mixed,'v5_auxiliary':d.group.eq('RNA_candidate'),'v5_primary':d.group.eq('Treg')}.items():
  rows.append((s,k,v.sum()))
 d.loc[f|d.group.eq('RNA_candidate')|d.id.eq(106665 if s=='Cervilla' else 49654),['id','x','y','mixed','group','pTreg','T_width','FOXP3']].to_parquet(R/s/'priority_old.parquet',index=False)
 tr=js(N/s/'metadata/transform.json');sc=np.array(tr['local_to_model_scale_xy']);origin=np.array(tr['fullres_crop_xyxy'][:2]);pick=[('small_nucleus',g & d.nuclear_area_um2.lt(25)),('epithelium',d.epi_eligible),('interface',g & d.mixed),('uncertain',d.geometry_review & d.in_D)]
 for i,(k,b) in enumerate(tr['representative_model_bboxes'].items()):fields.append((s,k+'_prior',(b[0]+b[2])//2,(b[1]+b[3])//2))
 for k,sel in pick:
  for z in d[sel].sample(n=2,random_state=123).itertuples():
   x,y=np.rint((np.array([z.x,z.y])/tr['native_mpp']-origin+.5)*sc-.5).astype(int);fields.append((s,k+'_'+str(z.id),x,y))
tab(pd.DataFrame(rows,columns=['section','condition','n']),'v5_funnel.tsv')
tab(pd.DataFrame(fields,columns=['section','field','mx','my']),'calibration_fields.tsv')
m=rd(V5/'input_manifest.tsv');m['readable']=m.path.map(lambda p:Path(p).exists());put(m,R/'input_manifest.tsv')
