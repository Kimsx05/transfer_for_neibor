from base6 import *
notes=rd(R/'tables/visual_notes.tsv');assert not notes.duplicated(['section','id']).any()
for s in sections:
 o=R/s;d=pd.read_parquet(o/'identity_frozen.parquet');n=notes[notes.section==s].set_index('id')
 d['review_action']=d.id.map(n.action).fillna('pending');bad=d.review_action.isin(['assignment_uncertain','unusable']);d.loc[bad,'mixed_v6']='assignment_uncertain'
 for col in ['A','B','main','strict','extended','control','epi']:d.loc[bad,col]=False
 d.loc[bad,'geometry_assignable']=False;d.loc[bad,'outcome']='assignment_uncertain'
 d.to_parquet(o/'cells_pre_MP10.parquet',index=False)
 q=pd.read_parquet(o/'review_queue.parquet');q['image']=q.review_id.map(lambda k:'cards/'+k.replace(':','_')+'.jpg');q['rendered']=q.image.isin(['cards/'+p.name for p in (o/'cards').glob('*.jpg')]);q['visual_note']=q.id.map(n.observation);q['review_state']=np.where(q.visual_note.notna(),'opened_reviewed','pending')
 aux=R/'tables/unlocalized_visual.tsv'
 if aux.exists():
  t=rd(aux);f=t[t.section.eq(s)].set_index('review_id').observation
  mapped=q.review_id.map(f);use=mapped.notna();q.loc[use,'visual_note']=mapped[use];q.loc[use,'review_state']='opened_reviewed'
 q.to_parquet(o/'review_index.parquet',index=False)
 tr=rd(o/'priority_trace.tsv');tr['final_after_review']=tr.new_id.map(d.set_index('id').outcome);put(tr,o/'priority_final.tsv')
notes['file']=notes.apply(lambda z:str(R/z.section/'cards'/('cell_'+str(z.id)+'.jpg')),axis=1);notes['review_stage']='identity_MP10_hidden';notes['reviewer']='Codex AI';put(notes,R/'tables/visual_review.tsv')
