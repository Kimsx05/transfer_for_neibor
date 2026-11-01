"""One seeded, local binary min-cut on observed Visium edges."""
from pathlib import Path
import json
import sys
from collections import Counter

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/mp10_pathology_regions_v3_3'
sys.path.insert(0,str(BASE/'scripts'))
import run_v3_3 as geo
import register_author as paper
from run_v3_3 import np, pd, M, graph, components, Image
from scipy.ndimage import label, uniform_filter
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow, breadth_first_order

CONFIG = json.loads((OUT/'config.json').read_text())
INSIDE, OUTSIDE = 1, -1

def load_section(sec):
    return pd.read_parquet(BASE/'by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_3.parquet')

def registered_parents(sec, d, adj):
    """Match un-eroded printed parent components to fixed cores on the observed graph."""
    folder = BASE/'registration'/sec
    meta = json.loads((folder/'registration.json').read_text())
    if meta['status'] != 'reliable_coarse_HE_registration':
        return None
    refs = pd.read_csv(folder/'coarse_reference.tsv.gz',sep='\t').set_index('barcode').loc[d.barcode]
    assert np.array_equal(refs.author_coarse_reference,d.author_coarse_reference)
    x,y,w,h = meta['annotation_crop']
    whole = np.asarray(Image.open(paper.PRE/'sources/image1.png').convert('RGB'))
    crop = whole[y:y+h,x:x+w]
    sep = np.linalg.norm(crop[:,:,None,:]-paper.palette(whole)[None,None,:,:],axis=3)
    colors = sep.argmin(2)
    exact = sep.min(2)<28
    coords = np.rint(refs[['author_x','author_y']].to_numpy()).astype(int)
    xx,yy = coords.T
    on = (xx>=0)&(xx<w)&(yy>=0)&(yy<h)
    def sample(img):
        vals = np.zeros(len(d),dtype=img.dtype)
        vals[on]=img[yy[on],xx[on]]
        return vals
    core = np.select([refs.author_coarse_reference=='author_coarse_inside',
                      refs.author_coarse_reference=='author_coarse_outside'],[INSIDE,OUTSIDE],default=0)
    strength = {side:np.zeros(len(d)) for side in [INSIDE,OUTSIDE]}
    parents = {side:np.zeros(len(d),int) for side in [INSIDE,OUTSIDE]}
    centers = {}
    parent_image = np.zeros((h,w),int)
    win = 2*int(np.ceil(meta['exclusion_buffer_author_px']))+1
    known = d.lesion_location.isin(['inside','outside']).to_numpy()
    for side,classes in [(INSIDE,[0,1]),(OUTSIDE,[6,7])]:
        raster, _ = label(exact & np.isin(colors,classes))
        center = sample(raster)
        centers[side] = center
        for parent_id in np.unique(center[core==side]):
            assert parent_id>0
            mask = raster==parent_id
            parent_image[mask]=side
            frac = sample(uniform_filter(mask.astype(float),size=win,mode='constant'))
            seeds = (core==side)&(center==parent_id)
            eligible = known & ((frac>=CONFIG['minimum_parent_neighborhood_fraction'])|seeds) & (core!=-side)
            # A raster component can intersect disjoint pieces of observed tissue.
            # Only the connected part containing a matching hard seed is eligible.
            lab, _ = components(adj,eligible)
            reachable = eligible & np.isin(lab,np.unique(lab[seeds]))
            prefer = reachable & (frac>strength[side])
            strength[side][prefer]=frac[prefer]
            parents[side][prefer]=parent_id
    return core, strength, parents, centers, crop, parent_image

def supported_compartments(d, adj, loc, ref):
    core, strength, parents, centers = ref[:4]
    _, groups = components(adj,loc=='outside')
    results=[]
    for ids in groups:
        if not all(len(adj[i])==6 for i in ids):
            continue
        rim = {j for i in ids for j in adj[i] if loc[j]!='outside'}
        if not rim or any(loc[j]!='inside' for j in rim):
            continue
        ring = np.array(sorted(rim))
        if not (d.v3_1_lesion_location.iloc[ring].eq('outside').to_numpy() &
                ((core[ring]==INSIDE)|(strength[INSIDE][ring]>=0.5))).any():
            continue
        membership=centers[INSIDE][ids]
        eligible=(membership>0)&(membership==parents[INSIDE][ids])
        counts=Counter(membership[eligible])
        if not counts:
            continue
        parent, count=counts.most_common(1)[0]
        if (count>len(ids)*CONFIG['enclosure_single_parent_fraction_exclusive'] and (parents[INSIDE][ring]==parent).any()
            and not (core[ids]==OUTSIDE).any() and not (strength[OUTSIDE][ids]>=0.5).any()
            and not (centers[OUTSIDE][ids]>0).any()):
            results.append((ids,int(parent)))
    return results

def mincut(d, adj, base, ref, enclosed):
    core,strength,parents=ref[:3]
    n=len(d)
    frozen=~d.lesion_location.isin(['inside','outside']).to_numpy()
    # Unary costs. Prior and molecular terms are bounded and symmetric.
    prior=np.where(base=='inside',1.,-1.)
    local=d.domain_local_evidence_score.to_numpy(float)
    molecular=CONFIG['molecular_weight']*np.clip(np.nan_to_num(local-.42),-1,1)
    preference=CONFIG['prior_weight']*prior+molecular+CONFIG['parent_weight']*(strength[INSIDE]-strength[OUTSIDE])
    active=(strength[INSIDE]>0)|(strength[OUTSIDE]>0)
    halo={j for i in np.flatnonzero(active) for j in adj[i]}
    active[list(halo)]=True
    active &= ~frozen
    hard=np.zeros(n,int)
    hard[~active]=np.where(base[~active]=='inside',INSIDE,OUTSIDE)
    hard[enclosed]=INSIDE
    hard[core!=0]=core[core!=0]
    # Source-side = inside. Capacities are integral for scipy's max-flow.
    src,sink=n,n+1
    row=[]; col=[]; cap=[]
    def edge(i,j,w):
        if w>0:
            row.append(i);col.append(j);cap.append(int(round(w*1000)))
    dom=d.domain_component_id.to_numpy()
    for i in range(n):
        edge(src,i,max(preference[i],0)+(CONFIG['hard_constraint_cost'] if hard[i]==INSIDE else 0))
        edge(i,sink,max(-preference[i],0)+(CONFIG['hard_constraint_cost'] if hard[i]==OUTSIDE else 0))
        if frozen[i]:
            continue
        for j in adj[i]:
            if j>i and not frozen[j]:
                weight=CONFIG['within_domain_edge_cost'] if dom[j]==dom[i] else CONFIG['between_domain_edge_cost']
                edge(i,j,weight);edge(j,i,weight)
    capacities=coo_matrix((cap,(row,col)),shape=(n+2,n+2),dtype=np.int32).tocsr()
    result=maximum_flow(capacities,src,sink)
    residual=(capacities-result.flow).tocsr()
    residual.eliminate_zeros()
    reach=breadth_first_order(residual,src,directed=True,return_predecessors=False)
    loc=base.copy()
    loc[active]='outside'
    observed = reach[reach < n]
    loc[observed[active[observed]]]='inside'
    assert all(loc[core==INSIDE]=='inside') and all(loc[core==OUTSIDE]=='outside')
    return loc,active

def segment(sec,d):
    base=d.lesion_location.to_numpy().copy()
    adj,_=graph(d)
    ref=registered_parents(sec,d,adj) if sec.startswith('GSE319536_') else None
    reason=np.full(len(d),'carry_v3.3_candidate',object)
    support=d.position_support.to_numpy().copy()
    parent=np.zeros(len(d),int)
    if ref is None:
        return base,reason,support,parent,None,[]
    core,strength,parents=ref[:3]
    initial=supported_compartments(d,adj,base,ref)
    enclosed=np.concatenate([ix for ix,p in initial]).astype(int) if initial else np.array([],int)
    loc,active=mincut(d,adj,base,ref,enclosed)
    # Reject newly detached islands with no same-side prior position or registered seed.
    # This connectivity condition does not depend on island size.
    proposed = loc.copy()
    for side, target in [(INSIDE, 'inside'), (OUTSIDE, 'outside')]:
        _, regions = components(adj, proposed == target)
        for ids in regions:
            if not (base[ids] == target).any() and not (core[ids] == side).any():
                loc[ids] = base[ids]
    # Conditional region assignment, not a mask-wide hole fill.
    later=supported_compartments(d,adj,loc,ref)
    for ix,p in later:
        loc[ix]='inside'
    changed=loc!=base
    reason[changed]='seeded_parent_and_observed_domain_continuity'
    support[changed]='author_parent_spatial_inference'
    for side,target in [(INSIDE,'inside'),(OUTSIDE,'outside')]:
        selected=active&(loc==target)&(parents[side]>0)
        parent[selected]=parents[side][selected]*side
    reason[core!=0]='retain_registered_core_constraint'
    for ix,p in initial+later:
        reason[ix]='enclosed_observed_compartment_same_tumour_parent'
        support[ix]='author_parent_and_complete_enclosure'
        parent[ix]=p
    return loc,reason,support,parent,ref,initial+later
