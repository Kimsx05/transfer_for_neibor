"""Attach the final decision and the actual viewing state to each source-derived review task."""
from base6 import *
for section in sections:
 o=R/section;q=pd.read_parquet(o/'review_index.parquet');d=pd.read_parquet(o/'cells_scored.parquet').set_index('id',verify_integrity=True);src=pd.read_parquet(o/'sources.parquet');links=pd.read_parquet(o/'source_links.parquet')
 q['outcome_before_review']=q.outcome;q['final_outcome']=q.id.map(d.outcome).fillna('unlocalized')
 q['pending_visual']=q.review_state.eq('pending');q['archive_path']=q.archive_shard.map(lambda i:f'card_archive/{i:04d}.sqlite')
 q['evidence_status']=np.select([q.pending_visual,q.final_outcome.isin(['unlocalized','assignment_uncertain','low_UMI','insufficient_subtype','no_direct_T','weak_major_T'])],['awaiting_visual_review','opened_identity_unconfirmed'],default='opened_with_disposition')
 q.to_parquet(o/'review_index.parquet',index=False)
 cv=q[q.id.notna()].set_index('id')
 links['final_outcome']=links.id.map(d.outcome);links['card_seen']=links.id.map(cv.review_state).eq('opened_reviewed')
 assert links.final_outcome.notna().all()
 f=links.groupby('unit').agg(linked_cells=('id','size'),linked_cards_opened=('card_seen','sum'))
 op=links.groupby(['unit','final_outcome']).size().unstack(fill_value=0).astype('int32').add_prefix('n_')
 summary=src[['unit','n_final_cells','roi']].set_index('unit').join(f).join(op).fillna(0)
 summary['unlocalized_card_seen']=q.set_index('review_id').review_state.reindex(summary.index).eq('opened_reviewed') & summary.n_final_cells.eq(0)
 assert (summary.linked_cells==summary.n_final_cells).all()
 summary.reset_index().to_parquet(o/'source_dispositions.parquet',index=False)
