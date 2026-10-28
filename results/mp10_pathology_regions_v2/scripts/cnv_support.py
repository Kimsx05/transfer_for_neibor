import os,sys,json,time
from pathlib import Path
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
R=Path(__file__).resolve().parents[2];O=R/'mp10_pathology_regions_v2'
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d
m=pd.read_csv(O/'tables/section_manifest.tsv',sep='\t');pilot=json.loads((O/'tables/pilot_sections.json').read_text());order=pilot+[s for s in m.section_id if s not in pilot]
def pearson(a,b):
 a=a-a.mean(axis=0);b=b-b.mean();den=np.sqrt((a*a).sum(0)*(b*b).sum());return np.divide((a*b[:,None]).sum(0),den,out=np.zeros(a.shape[1]),where=den>0)
for sec in order:
 r=m.set_index('section_id').loc[sec];dest=Path(r.output_dir)
 if (dest/'cnv_v3.done').exists():continue
 print(time.strftime('%H:%M:%S'),sec,flush=True)
 try:
  meta=pd.read_csv(dest/'input_spots.tsv',sep='\t');ar=np.load(r.cna_path,allow_pickle=False);A=ar['profiles'];bc=ar['barcodes'].astype(str);chrom=ar['chrom'];pos=ar['chrompos'];assert len(bc)==A.shape[1] and len(chrom)==A.shape[0];meta=meta.set_index('barcode');mo=meta.loc[bc];cl=pd.read_csv(r.cluster_path,sep='\t').set_index('barcode').loc[bc];clusters=cl.yan_adapter_cluster.to_numpy();ids=sorted(set(clusters));cons=np.stack([np.median(A[:,clusters==c],axis=1) for c in ids],axis=1);globalcons=cons.mean(1);cormean=pearson(cons,globalcons);msq=(A*A).mean(0);cmsq=(cons*cons).mean(0)
  # Windows preserve native units; only a separate reference-relative contrast is formed.
  window_keys=sorted(set(zip(chrom[chrom<=22],(pos[chrom<=22]//20000000))));windowbins=[np.flatnonzero((chrom==c)&(pos//20000000==w)) for c,w in window_keys];use=[i for i,x in enumerate(windowbins) if len(x)>=10 and np.ptp(pos[x])>=10000000];windowbins=[windowbins[i] for i in use];window_keys=[window_keys[i] for i in use];wc=np.array([np.median(cons[ix],axis=0) for ix in windowbins]);wp=np.array([np.median(A[ix],axis=0) for ix in windowbins]);baseline=(mo.native_copykat_label=='diploid').to_numpy();refreason='native high-confidence diploid spots'
  if baseline.sum()<50:
   baseline=msq<=np.quantile(msq,.25);refreason='lowest-burden quartile technical reference; not certified normal'
  ref=np.median(wp[:,baseline],axis=1);noise=1.4826*np.median(np.abs(wp[:,baseline]-ref[:,None]),axis=1);events=[];cr=[];sp=[]
  for j,c in enumerate(ids):
   ix=clusters==c;n=int(ix.sum());cor=pearson(A[:,ix],cons[:,j]);coherence=float(np.median(cor));threshold=np.maximum(.04,3*noise/np.sqrt(max(1,n)));delta=wc[:,j]-ref;global_offset=float(np.median(delta));structured_delta=delta-global_offset;ev=(np.abs(structured_delta)>=threshold);evchr=sorted(set(window_keys[k][0] for k in np.flatnonzero(ev)));native=mo.loc[ix,'native_copykat_label'];an=float(native.eq('aneuploid').mean());low=float(native.eq('c1:diploid:low.conf').mean());quality=mo.baseline_confidence.iloc[0];baseline_uncertain='UNCERTAIN' in str(quality)
   # Strong requires coherent multichromosome departures; native label is additional, not sole evidence.
   if n>=10 and coherence>=.3 and (len(evchr)>=2 or (len(evchr)>=1 and an>=.5)):support='strong'
   elif n>=5 and coherence>=.2 and len(evchr)>=1:support='moderate'
   else:support='weak_or_no_detectable_pattern'
   if support=='strong' and an<.5:support='coherent_pattern_native_diploid_conflict'
   if baseline_uncertain and support=='strong':support='moderate_baseline_caution'
   for k in np.flatnonzero(ev):
    bins=windowbins[k];events.append(dict(cluster=c,chromosome=window_keys[k][0],start=int(pos[bins].min()),end=int(pos[bins].max()),n_bins=len(bins),native_consensus_median=float(wc[k,j]),reference_median=float(ref[k]),reference_relative_delta=float(delta[k]),threshold_native_units=float(threshold[k]),direction='relative_gain_after_removing_global_offset' if structured_delta[k]>0 else 'relative_loss_after_removing_global_offset',structured_delta=float(structured_delta[k]),global_offset=global_offset))
   cr.append(dict(cluster=c,n_spots=n,global_offset_relative_to_reference=global_offset,chromosome_specific_range=float(np.ptp([np.median(delta[[kk for kk,key in enumerate(window_keys) if key[0]==cc]]) for cc in sorted(set(key[0] for key in window_keys))])),consensus_MSQ=float(cmsq[j]),mean_spot_MSQ=float(msq[ix].mean()),consensus_correlation_to_unweighted_cluster_mean=float(cormean[j]),median_spot_consensus_correlation=coherence,n_event_chromosomes=len(evchr),event_chromosomes=','.join(map(str,evchr)),n_event_windows=int(ev.sum()),native_aneuploid_fraction=an,native_low_conf_diploid_fraction=low,CNV_support=support,baseline_confidence=quality,reference_description=refreason,n_reference_spots=int(baseline.sum()),rule='coherent >=10Mb window deviations; abs(delta-median(delta))>=max(0.04,3*reference_MAD/sqrt(cluster_n)); no Yan MSQ threshold transfer'))
   for k,ccor in zip(np.flatnonzero(ix),cor):sp.append(dict(barcode=bc[k],CNA_cluster=c,CNV_support=support,n_event_chromosomes=len(evchr),native_spot_MSQ=float(msq[k]),cluster_consensus_MSQ=float(cmsq[j]),spot_consensus_correlation=float(ccor),baseline_confidence=quality))
  pd.DataFrame(cr).to_csv(dest/'cnv_cluster_support.tsv',sep='\t',index=False);pd.DataFrame(sp).to_csv(dest/'cnv_spot_support.tsv',sep='\t',index=False);pd.DataFrame(events,columns=['cluster','chromosome','start','end','n_bins','native_consensus_median','reference_median','reference_relative_delta','threshold_native_units','direction','structured_delta','global_offset']).to_csv(dest/'cnv_chromosome_events.tsv',sep='\t',index=False)
  # Compact true profile display: order spots by original CNA cluster; deterministically sample <=600.
  ordered=np.argsort(clusters,kind='stable');chosen=ordered[np.linspace(0,len(ordered)-1,min(600,len(ordered))).astype(int)];bins=np.arange(0,len(chrom),max(1,int(np.ceil(len(chrom)/1800))));np.savez_compressed(dest/'cnv_plot_profiles.npz',profiles=A[bins][:,chosen].astype('float32'),chrom=chrom[bins],chrompos=pos[bins],spot_clusters=clusters[chosen],consensus=cons[bins].astype('float32'),cluster_ids=ids,barcodes=bc[chosen]);(dest/'cnv_v3.done').write_text('complete; actual native profile read')
  del A,ar,wp,cons
 except Exception as e:
  import traceback;(dest/'cnv.error').write_text(traceback.format_exc());print('ERROR',sec,str(e),flush=True)
