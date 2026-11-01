"""Build one self-contained section delivery from the cached v3.3 inputs."""
import gzip
import tempfile
from concurrent.futures import ProcessPoolExecutor
from regions import *
from run_v3_3 import exterior_empty, FULL, OFF, VERT

BAND_NAMES = {
    'Interior (>w confirmed)': 'candidate_inside_beyond_w',
    'Outside (>w confirmed)': 'candidate_outside_beyond_w',
    'Inner margin': 'inner_margin', 'Outer margin': 'outer_margin',
    'Band unresolved': 'band_unresolved', 'Location unresolved': 'unlocalized',
    'QC excluded': 'quality_excluded',
    'Inside; outer edge unobserved': 'inside_outer_unobserved',
    'Outside; no usable lesion edge': 'outside_outer_unobserved',
}
DIST_FIELDS = ['signed_distance_um','abs_distance_lower_bound_um','abs_distance_upper_bound_um',
               'distance_to_observed_outer_candidate_um','distance_to_unobserved_edge_um']
WIDTHS = (100,200,300)

def section_dir(sec):
    return OUT/'by_dataset'/M.loc[sec,'dataset_id']/sec

def final_ids(values):
    return values.astype(str).str.replace(':v3_3:',':F:',regex=False)

def handle_internal_outside(d, sec, adj, lookup, edges):
    """Only an externally accessible outside group receives macro-outer distance."""
    loc=d.lesion_location.to_numpy()
    ext=exterior_empty(lookup)
    full=FULL[sec]
    rc=list(zip(d.array_row,d.array_col))
    d['region_role']=np.where(loc=='inside','candidate_lesion_body',np.where(loc=='outside','candidate_external_tissue',loc))
    enclosed=d.final_reason.eq('enclosed_observed_compartment_same_tumour_parent')
    d.loc[enclosed,'region_role']='lesion_internal_compartment'
    _,groups=components(adj,loc=='outside')
    inaccessible=[]
    for ix in groups:
        accessible=any((q:=(rc[i][0]+a,rc[i][1]+b)) in ext and q not in full for i in ix for a,b in OFF)
        if accessible:
            continue
        rim={j for i in ix for j in adj[i] if loc[j]!='outside'}
        closed=all(len(adj[i])==6 for i in ix) and bool(rim) and all(loc[j]=='inside' for j in rim)
        d.loc[ix,'region_role']='enclosed_outside_compartment' if closed else 'outside_external_route_unobserved'
        inaccessible.extend(ix)
    for w in WIDTHS:
        col=f'distance_band_{w}um'
        mapped=d[col].map(BAND_NAMES)
        assert mapped.notna().all()
        d[col]=mapped
        d.loc[inaccessible,col]='outside_no_external_route'
        d.loc[inaccessible,f'band_status_{w}um']='no_observed_external_route'
    d.loc[inaccessible,DIST_FIELDS]=np.nan
    d.loc[inaccessible,['reference_lesion_roi_id','reference_boundary_segment_id']]=''
    d.loc[inaccessible,'exact_distance_status']='no_observed_external_route'
    d.loc[inaccessible,'reference_boundary_support']='no_candidate_reference'
    outer=edges[edges.boundary_type=='candidate_observed_outer']
    assert not set(outer.outside_index)&set(inaccessible)
    d['observed_outer_contact']=d.index.isin(set(outer.inside_index)|set(outer.outside_index))
    inner=edges[edges.boundary_type=='internal_nonlesion_compartment_contact']
    d['internal_compartment_contact']=d.index.isin(set(inner.inside_index)|set(inner.outside_index))
    return d

def split_table(df):
    # Portable column/row table; no object pickle or upstream expression matrix.
    return json.loads(df.to_json(orient='split',index=False,double_precision=15))

def save_geometry(src, out, d, e, un, co, sec):
    rois=json.loads((src/'lesion_rois_um.geojson').read_text())
    for f in rois['features']:
        f['properties']['lesion_roi_id']=f['properties']['lesion_roi_id'].replace(':v3_3:',':F:')
    e=e.copy()
    e['inside_barcode']=d.barcode.to_numpy()[e.inside_index.astype(int)]
    e['outside_barcode']=[d.barcode.iloc[j] if j>=0 else '' for j in e.outside_index]
    e['segment_id']=final_ids(e.segment_id)
    e['lesion_roi_id']=final_ids(e.lesion_roi_id)
    gaps=pd.read_csv(src/'internal_gap_assessment.tsv',sep='\t')
    old=BASE/'by_dataset'/M.loc[sec,'dataset_id']/sec
    nests=pd.read_csv(old/'tumor_nest_contact_edges.tsv.gz',sep='\t')
    nests['inside_lesion_roi']=[bool(d.lesion_location.iloc[i]=='inside' and d.lesion_location.iloc[j]=='inside') for i,j in zip(nests.inside_index,nests.outside_index)]
    bundle={'section_id':sec,'coordinate_scale':co,'roi_polygons_um':rois,
            'roi_edges':split_table(e),'unknown_K_edges':split_table(un),
            'unobserved_gap_cells':split_table(gaps),'secondary_nest_edges':split_table(nests)}
    with gzip.open(out/'geometry.json.gz','wt',encoding='utf8',compresslevel=5) as stream:
        json.dump(bundle,stream,separators=(',',':'),allow_nan=False)
    return e

def process(sec):
    out=section_dir(sec)
    out.mkdir(parents=True,exist_ok=True)
    old=BASE/'by_dataset'/M.loc[sec,'dataset_id']/sec
    d=load_section(sec)
    base=d.lesion_location.copy()
    loc,why,support,parent,ref,holes=segment(sec,d)
    changed=(loc!=base.to_numpy())
    co=json.loads((old/'coordinate_scale.json').read_text())
    aff=np.asarray(co['array_to_image_affine'])
    image=np.asarray(Image.open(M.loc[sec,'image_path']).convert('RGB'))
    adj,lookup=graph(d)
    tc,tcs=components(adj,np.ones(len(d),bool))
    d['v3_3_lesion_location']=base
    d['lesion_location']=loc
    d['final_reason']=why
    d['position_support']=support
    d['signed_author_parent_id']=parent
    with tempfile.TemporaryDirectory() as temp:
        if changed.any():
            src=Path(temp)
            d,e,un,poly,_=geo.new_geometry(d,adj,lookup,tc,tcs,aff,co,image,sec,src)
        else:
            src=old
            e=pd.read_csv(old/'lesion_boundary_edges.tsv.gz',sep='\t')
            un=pd.read_csv(old/'unknown_region_boundary.tsv.gz',sep='\t')
            keys=np.c_[3*d.array_col,3*d.array_row].astype(int)[:,None,:]+VERT
            uv=np.concatenate([keys[:,:,:1]/6,keys[:,:,1:]/(2*np.sqrt(3)),np.ones((len(d),6,1))],axis=2)
            poly=uv@aff
        d=handle_internal_outside(d,sec,adj,lookup,e)
        e=save_geometry(src,out,d,e,un,co,sec)
    for col in ['lesion_roi_id','reference_lesion_roi_id','reference_boundary_segment_id']:
        d[col]=final_ids(d[col])
    d['section_row_index']=np.arange(len(d))
    raw=d.author_panel_class_at_center
    raw_side=raw.map({'Tumor Cell':'inside','Tumor Stroma':'inside','Urothelium':'outside','Normal Stroma':'outside'})
    d['reference_state']=np.select([
        d.author_coarse_reference.isin(['author_coarse_inside','author_coarse_outside']),
        raw_side.notna() & raw_side.ne(d.lesion_location),
        (ref is not None)&raw_side.eq(d.lesion_location)],
        ['registered_core_constraint','possible_parent_reference_conflict','regional_context_agrees'],
        default='no_binary_reference')
    d['interface_status']='no_observed_candidate_outer_in_block'
    for k in e.loc[e.boundary_type.eq('candidate_observed_outer'),'tissue_component'].unique():
        d.loc[tcs[k],'interface_status']='observed_candidate_outer_in_block'
    fields=['dataset_id','section_id','barcode','x','y','array_row','array_col','x_um','y_um','section_row_index',
            'tissue_component_id','domain_component_id','v3_3_lesion_location','lesion_location','lesion_roi_id',
            'final_reason','position_support','author_panel_class_at_center','author_coarse_reference','reference_state',
            'signed_author_parent_id','region_role','interface_status','observed_outer_contact','internal_compartment_contact',
            'reference_lesion_roi_id','reference_boundary_segment_id','reference_boundary_support','exact_distance_status',
            *DIST_FIELDS,*[f'distance_band_{w}um' for w in WIDTHS],*[f'band_status_{w}um' for w in WIDTHS]]
    # Section labels are the source for streaming one complete compressed cohort table.
    d[fields].to_csv(out/'labels.tsv.gz',sep='\t',index=False)
    from figures import plot_section
    plot_section(sec,d,base,poly,image,e,co,out,ref,holes)
    return dict(dataset_id=M.loc[sec,'dataset_id'],section_id=sec,n_spots=len(d),position_changes=int(changed.sum()),
                inside=int((loc=='inside').sum()),outside=int((loc=='outside').sum()),
                unlocalized=int((loc=='unlocalized').sum()),candidate_edges=int(e.boundary_type.eq('candidate_observed_outer').sum()),
                precise_distance=int(d.signed_distance_um.notna().sum()),
                band_200_determined=int(d.distance_band_200um.isin(['candidate_inside_beyond_w','candidate_outside_beyond_w','inner_margin','outer_margin']).sum()),
                enclosed_assigned=int(d.region_role.eq('lesion_internal_compartment').sum()))

def run(sections):
    with ProcessPoolExecutor(max_workers=2) as pool:
        for result in pool.map(process,sections):
            print(json.dumps(result),flush=True)

if __name__=='__main__':
    run(sys.argv[1:] or list(M.index))
