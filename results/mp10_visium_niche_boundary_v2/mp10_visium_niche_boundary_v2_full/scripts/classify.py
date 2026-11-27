from lib import *
import argparse
p=argparse.ArgumentParser();p.add_argument('--config',default=str(O/'provisional_classification.json'));p.add_argument('--version',default='primary');args=p.parse_args()
if not args.version.replace('_','').replace('-','').isalnum():raise ValueError('version')
out=O/'classification' if args.version=='primary' else O/'versions'/args.version/'classification'
out.mkdir(exist_ok=True,parents=True);cfg=json.loads(Path(args.config).read_text());signature=sha(Path(args.config))
lock=out/'classification_lock.json'
if lock.exists() and json.loads(lock.read_text())['config_sha256']!=signature:raise RuntimeError('Use a new --version for new thresholds')
w=read(O/'classification/immune_windows.tsv.gz');links=read(O/'classification/window_links.tsv.gz');nat=read(O/'classification/natural_clusters.tsv').set_index('section_id');x=read(O/'classification/immune_input.tsv.gz').set_index('section_id')
assert not any('mp10' in col.lower() or 'p95' in col.lower() for col in w)
rows=[];win=[]
for (sec,level),v in w.groupby(['section_id','level']):
    v=v.copy();width=int(level[6:]);minE=cfg['min_E_'+str(width)];minout=cfg['min_outer_'+str(width)]
    ok=v.E_n>=minE;per=v.per_n>=minout;density=v.E_T_density;frac=v.E_T_fraction;lag=v.T_neighbor_E
    # Exclusion requires sufficient absolute peripheral supply in addition to a low ratio.
    excluded=ok&per&(v.out_T_density>=cfg['T_floor'])&(v.P_T<cfg['entry_floor'])&(density<cfg['T_strong'])
    poor=ok&(density<cfg['T_floor'])&(frac<cfg['E_T_fraction_floor'])
    ins=ok&(density>=cfg['T_floor'])&(frac>=cfg['E_T_fraction_floor'])&(lag>=cfg['neighbor_floor'])&((per&(v.P_T>=cfg['entry_floor']))|(~per&(density>=cfg['T_strong'])&(lag>=cfg['neighbor_strong'])))
    v['state_reason']=np.select([~ok,excluded,poor,ins],['low_E_support','excluded-like','immune-poor','local_entry_or_internal_mixing'],default='indeterminate')
    v['local_state']=np.select([ins,excluded|poor],['infiltrated','non_infiltrated'],default='unresolved')
    ll=links[(links.section_id==sec)&(links.level==level)];states=v.set_index('local_id').local_state
    matched=ll[ll.a.map(states)==ll.b.map(states)]
    v['neighbor_support']=v.local_id.isin(set(matched.a)|set(matched.b)) & (v.local_state!='unresolved')
    win.append(v);total=v.E_n.sum();usable=v.loc[ok,'E_n'].sum()
    masks=[(v.local_state==s)&v.neighbor_support for s in ['infiltrated','non_infiltrated']]
    n1,n0=[int(a.sum()) for a in masks];a1,a0=[v.loc[a,'E_n'].sum() for a in masks]
    f1=safe(a1,usable);f0=safe(a0,usable);minw=cfg['minimum_state_windows_'+str(width)]
    covered=safe(a1+a0,total);sh=np.array([a1,a0],float);probs=sh[sh>0]/sh.sum() if sh.sum() else []
    for rule,threshold in [('primary',cfg['mixed_area_share']),('loose',cfg['loose_share']),('strict',cfg['strict_share'])]:
        if not usable:group='not_evaluable';reason='no_supported_E_windows'
        elif covered<cfg['min_E_area_coverage']:group='uncertain';reason='insufficient_spatially_supported_area'
        elif min(n1,n0)>=minw and min(f1,f0)>=threshold:group=G[2];reason='two_coherent_local_patterns'
        elif n1>=minw and f1>=cfg['dominant_area_share']:group=G[0];reason='supported_entry_predominates'
        elif n0>=minw and f0>=cfg['dominant_area_share']:
            group=G[1];z=v[masks[1]].groupby('state_reason').E_n.sum();reason=z.idxmax()
        else:group='uncertain';reason='no_dominant_or_dual_supported_pattern'
        rows.append(dict(section_id=sec,level=level,rule=rule,provisional_architecture_group=group,group_reason=reason,confidence='limited' if group not in G else 'higher_coverage' if covered>=.7 else 'moderate_coverage',evaluability='classified' if group in G else group,mode='normal_projection' if M.at[sec,'source_context']=='adjacent_normal_source' else 'tumor_discovery',n_entry=n1,n_non=n0,E_area_eligible=usable*AREA,F_entry=f1,F_non=f0,F_unresolved=1-f1-f0,window_share_entry=safe(n1,sum(ok)),window_share_non=safe(n0,sum(ok)),effective_area_coverage=covered,H_supported=-sum(val*np.log2(val) for val in probs) if len(probs) else np.nan,**nat.loc[sec].to_dict()))
full=pd.DataFrame(rows);save(full,out/'section_classification_all.tsv');primary=full.query("rule=='primary' and level=='window500'").copy()
for label,q in [('scale1000',"rule=='primary' and level=='window1000'"),('loose',"rule=='loose' and level=='window500'"),('strict',"rule=='strict' and level=='window500'")]:
    primary[label]=primary.section_id.map(full.query(q).set_index('section_id').provisional_architecture_group)
save(primary,out/'section_classification.tsv');save(pd.concat(win),out/'window_states.tsv.gz')
lock.write_text(json.dumps(dict(config_sha256=signature,config=cfg,immune_table_hash=sha(O/'classification/immune_input.tsv.gz'),window_table_hash=sha(O/'classification/immune_windows.tsv.gz'),class_sha256=sha(out/'section_classification.tsv'),natural_model=json.loads((O/'classification/natural_model.json').read_text())),indent=2))
print(primary.groupby(['mode','provisional_architecture_group']).size())

