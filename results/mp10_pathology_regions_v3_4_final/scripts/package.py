"""Combine the cohort once, verify geometry against member keys, then publish local indexes."""
import gzip
import html
import json
import os
from html.parser import HTMLParser
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.parse import unquote

import shapely
from finalize import *
from figures import plt, LOCATION, save

MASTER=OUT/'final_spot_labels.tsv.gz'
NEAR=['inner_margin','outer_margin']
FAR=['candidate_inside_beyond_w','candidate_outside_beyond_w']
STYLE='<style>body{font:16px system-ui,sans-serif;color:#27303a;max-width:1240px;margin:30px auto;padding:0 20px}img{width:100%;height:auto}a{color:#164d75}td,th{padding:.4em .6em;border-bottom:1px solid #ddd;text-align:left}th{background:#f4f5f7}p{line-height:1.55}</style>'
def page(title,body):
    return f'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{html.escape(title)}</title>{STYLE}<body><h1>{html.escape(title)}</h1>{body}</body></html>'

def collect():
    cross=pd.read_csv(BASE/'tables/preserved_identity_crosswalk.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('section_id')
    existing=None
    frames=[]
    for sec in M.index:
        frag=section_dir(sec)/'labels.tsv.gz'
        if frag.exists():
            d=pd.read_csv(frag,sep='\t',float_precision='round_trip')
        else:
            if existing is None:
                existing=pd.read_csv(MASTER,sep='\t',float_precision='round_trip')
            d=existing[existing.section_id.eq(sec)].copy()
        d['analysis_unit_uid']=cross.at[sec,'analysis_unit_uid_v3_1']
        d['patient_uid']=cross.at[sec,'patient_uid']
        d['distance_unit']='um'
        d['region_version']='v3.4_final'
        if sec=='GSE238145_B16':
            d.loc[d.barcode.eq('AGATAACACGGTATAA-1_4'),'final_reason']='retain_v3.3_outside_six_observed_outside_neighbors'
        frames.append(d)
    data=pd.concat(frames,ignore_index=True)
    assert len(data)==265002 and data.section_id.nunique()==63 and data.dataset_id.nunique()==7
    assert not data.duplicated(['section_id','barcode']).any()
    temp=MASTER.with_name('final_spot_labels.tmp.gz')
    data.to_csv(temp,sep='\t',index=False)
    temp.replace(MASTER)
    return data

def as_frame(bundle,key):
    tab=bundle[key]
    return pd.DataFrame(tab['data'],columns=tab['columns'])

def validate(sec,d):
    b=load_section(sec)
    d=d.reset_index(drop=True)
    assert np.array_equal(b.barcode,d.barcode)
    for col in ['x','y','array_row','array_col','x_um','y_um','dataset_id','section_id','tissue_component_id']:
        assert np.array_equal(b[col],d[col]),(sec,col)
    assert np.array_equal(b.lesion_location,d.v3_3_lesion_location)
    changed=d.lesion_location.ne(b.lesion_location)
    if not sec.startswith('GSE319536_') or sec.endswith(('BC05','BC17')):
        assert not changed.any()
    excluded=b.lesion_location.isin(['unlocalized','quality_excluded'])
    assert d.loc[excluded,'lesion_location'].equals(b.loc[excluded,'lesion_location'])
    with gzip.open(section_dir(sec)/'geometry.json.gz','rt') as stream:
        bundle=json.load(stream)
    e=as_frame(bundle,'roi_edges')
    adj,lookup=graph(d)
    lab,groups=components(adj,d.lesion_location.eq('inside').to_numpy())
    polygons=bundle['roi_polygons_um']['features']
    assert len(polygons)==len(groups)==d.lesion_roi_id.dropna().nunique()
    members={rid:z.index for rid,z in d[d.lesion_location.eq('inside')].groupby('lesion_roi_id')}
    for obj in polygons:
        rid=obj['properties']['lesion_roi_id']
        ix=members[rid]
        assert len(ix)==obj['properties']['n_spots'] and len(set(lab[ix]))==1
        geom=shapely.geometry.shape(obj['geometry'])
        shapely.prepare(geom)
        assert shapely.covers(geom,shapely.points(d.x_um.iloc[ix],d.y_um.iloc[ix])).all()
    lookup_segment=e.set_index('segment_id')
    outer=e[e.boundary_type.eq('candidate_observed_outer')]
    for q in e.itertuples():
        i=int(q.inside_index); j=int(q.outside_index)
        assert d.barcode.iloc[i]==q.inside_barcode and d.lesion_roi_id.iloc[i]==q.lesion_roi_id
        if j>=0:
            assert j in adj[i] and d.tissue_component_id.iloc[j]==d.tissue_component_id.iloc[i]
            assert d.barcode.iloc[j]==q.outside_barcode
        if q.boundary_type=='candidate_observed_outer':
            assert j>=0 and d.lesion_location.iloc[j]=='outside'
    ordered=d[d.reference_boundary_segment_id.notna()]
    if len(ordered):
        matched=lookup_segment.loc[ordered.reference_boundary_segment_id]
        assert np.array_equal(ordered.reference_lesion_roi_id,matched.lesion_roi_id)
        assert matched.boundary_type.eq('candidate_observed_outer').all()
        xy=ordered[['x_um','y_um']].to_numpy()
        a=matched[['x1_um','y1_um']].to_numpy()
        ab=matched[['x2_um','y2_um']].to_numpy()-a
        t=np.clip(np.sum((xy-a)*ab,axis=1)/np.sum(ab*ab,axis=1),0,1)
        assert np.allclose(np.linalg.norm(xy-a-ab*t[:,None],axis=1),ordered.distance_to_observed_outer_candidate_um)
    exact=d.signed_distance_um.notna()
    assert exact.equals(d.exact_distance_status.eq('estimable'))
    assert d.loc[exact&d.lesion_location.eq('inside'),'signed_distance_um'].lt(0).all()
    assert d.loc[exact&d.lesion_location.eq('outside'),'signed_distance_um'].gt(0).all()
    for w in WIDTHS:
        band=d[f'distance_band_{w}um']
        near=band.isin(NEAR);far=band.isin(FAR)
        assert d.loc[near,'abs_distance_upper_bound_um'].le(w).all()
        assert d.loc[far,'abs_distance_lower_bound_um'].gt(w).all()
        assert d.loc[near|far,'reference_boundary_segment_id'].notna().all()
        assert d.loc[band.isin(['candidate_inside_beyond_w','inner_margin']),'lesion_location'].eq('inside').all()
        assert d.loc[band.isin(['candidate_outside_beyond_w','outer_margin']),'lesion_location'].eq('outside').all()
    core=d.author_coarse_reference
    assert d.loc[core.eq('author_coarse_inside'),'lesion_location'].eq('inside').all()
    assert d.loc[core.eq('author_coarse_outside'),'lesion_location'].eq('outside').all()
    # Every all-new ROI requires a prior or a same-side author seed. A seed-only
    # ROI cannot simply trace the buffered core within a larger observed block.
    if changed.any():
        for ix in groups:
            if b.lesion_location.iloc[ix].eq('outside').all():
                raise ValueError((sec,'new lesion disconnected from prior and core',ix[:4]))
            if b.v3_1_lesion_location.iloc[ix].eq('outside').all() and core.iloc[ix].eq('author_coarse_inside').all():
                block=d.tissue_component_id.iloc[ix[0]]
                assert len(ix)==int(d.tissue_component_id.eq(block).sum()),(sec,block,'buffer-core-only ROI')
    stats={'dataset_id':d.dataset_id.iloc[0],'section_id':sec,'spots':len(d),'changed':int(changed.sum()),
           'to_inside':int((changed&d.lesion_location.eq('inside')).sum()),'to_outside':int((changed&d.lesion_location.eq('outside')).sum()),
           'inside':int(d.lesion_location.eq('inside').sum()),'outside':int(d.lesion_location.eq('outside').sum()),
           'unlocalized':int(d.lesion_location.eq('unlocalized').sum()),'rois':len(groups),'candidate_edges':len(outer),
           'exact_distance':int(exact.sum()),'band_200_determined':int(d.distance_band_200um.isin(NEAR+FAR).sum()),
           'assigned_enclosed':int(d.region_role.eq('lesion_internal_compartment').sum()),
           'withdrawn_outside_route':int((b.reference_boundary_segment_id.fillna('').ne('') & d.exact_distance_status.eq('no_observed_external_route')).sum()),
           'old_margin':int(b.distance_band_200um.isin(['Inner margin','Outer margin']).sum()),
           'new_margin':int(d.distance_band_200um.isin(NEAR).sum()),
           'exited_margin':int((b.distance_band_200um.isin(['Inner margin','Outer margin'])&~d.distance_band_200um.isin(NEAR)).sum()),
           'entered_margin':int((~b.distance_band_200um.isin(['Inner margin','Outer margin'])&d.distance_band_200um.isin(NEAR)).sum())}
    cases={}
    if sec=='GSE319536_BC18':
        _,outside=components(adj,b.lesion_location.eq('outside').to_numpy())
        ix,= [ix for ix in outside if len(ix)==225 and (b.array_row.iloc[ix].min(),b.array_row.iloc[ix].max(),b.array_col.iloc[ix].min(),b.array_col.iloc[ix].max())==(14,36,23,57)]
        z=d.iloc[ix]
        assert all(len(adj[i])==6 for i in ix) and z.region_role.eq('lesion_internal_compartment').all()
        assert z.lesion_location.eq('inside').all() and z.distance_band_200um.eq('candidate_inside_beyond_w').all()
        assert b.distance_band_200um.iloc[ix].eq('Outside (>w confirmed)').sum()==224
        cases={'exact':int(z.signed_distance_um.notna().sum()),'classes':z.author_panel_class_at_center.value_counts().to_dict()}
    if sec=='GSE319536_BC13':
        ids=[rid for rid,z in b[b.lesion_location.eq('inside')].groupby('lesion_roi_id') if len(z) in [1,4,6,9] and z.v3_1_lesion_location.eq('outside').all()]
        core_rows=b.lesion_roi_id.isin(ids)
        oldmargin=b.reference_lesion_roi_id.isin(ids)&b.distance_band_200um.eq('Outer margin')
        assert len(ids)==5 and core_rows.sum()==21 and oldmargin.sum()==106
        assert set(['GATCCATAATCGCACG-1','GTACAGGACGTTACCT-1'])<=set(d.barcode[core_rows])
        assert d.loc[core_rows,'lesion_location'].eq('inside').all() and d.lesion_roi_id[core_rows].nunique()==1
        rid=d.lesion_roi_id[core_rows].iloc[0]
        assert d.distance_band_200um[oldmargin].eq('candidate_inside_beyond_w').all()
        cases={'merged_roi':rid,'merged_n':len(members[rid]),'former_margin':int(oldmargin.sum())}
    if sec=='GSE319536_BC14':
        single=(d.tissue_component_id.map(d.tissue_component_id.value_counts())==1)&b.lesion_location.eq('inside')&b.v3_1_lesion_location.eq('outside')
        assert single.any() and d.lesion_location[single].eq('inside').all() and d.signed_distance_um[single].isna().all()
        cases={'single':int(single.sum())}
    if sec=='GSE238145_B16':
        z=d[d.barcode.eq('AGATAACACGGTATAA-1_4')].iloc[0]
        assert z.lesion_location=='outside' and pd.isna(z.lesion_roi_id)
    return stats,cases

def write_indexes(data,summ):
    body='<p>63 张切片的正式计算候选位置。可用 section_id 与原始 barcode 一对一合并。</p><p><a href="final_spot_labels.tsv.gz">下载完整标签表（TSV.gz）</a> · <a href="METHODS.md">方法与读取示例</a> · <a href="RESULTS.md">结果说明</a></p>'
    (OUT/'overviews').mkdir(exist_ok=True)
    for dataset,frame in summ.groupby('dataset_id',sort=False):
        fig,axes=plt.subplots(int(np.ceil(len(frame)/4)),4,figsize=(14,3.4*int(np.ceil(len(frame)/4))),squeeze=False,layout='constrained')
        for ax in axes.flat:ax.axis('off')
        for ax,q in zip(axes.flat,frame.itertuples()):
            z=data[data.section_id.eq(q.section_id)]
            ax.scatter(z.x,z.y,s=2,c=[LOCATION[k][0] for k in z.lesion_location],linewidths=0)
            ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(f'{q.section_id}\n修改 {q.changed} 点 · 候选边段 {q.candidate_edges}',fontsize=8)
        fig.suptitle(dataset+'  |  红色为候选内侧，蓝色为候选外侧')
        save(fig,OUT/'overviews',dataset)
        body+=f'<h2>{dataset}</h2><a href="overviews/{dataset}.pdf">总览 PDF</a><img src="overviews/{dataset}.png" loading="lazy"><table><tr><th>切片</th><th>spot</th><th>修改点</th><th>候选外缘边段</th><th>200 μm 分带可判</th></tr>'
        for q in frame.itertuples():
            path=f'by_dataset/{dataset}/{q.section_id}'
            body+=f'<tr><td><a href="{path}/index.html">{q.section_id}</a></td><td>{q.spots}</td><td>{q.changed}</td><td>{q.candidate_edges}</td><td>{q.band_200_determined}</td></tr>'
            detail=f'<p>{q.spots} 个 spot。最终位置相对 v3.3 修改 {q.changed} 点。可观测候选外缘有 {q.candidate_edges} 条边段。</p><p><a href="../../../index.html">返回全队列</a> · <a href="../../../final_spot_labels.tsv.gz">完整标签表</a> · <a href="geometry.json.gz">本片几何数据</a></p>'
            for f in ['final','widths']+(['revision'] if q.changed else []):
                detail+=f'<p><a href="{f}.png">PNG</a> · <a href="{f}.pdf">PDF</a></p><img src="{f}.png" loading="lazy">'
            (OUT/path/'index.html').write_text(page(q.section_id,detail))
        body+='</table>'
    (OUT/'index.html').write_text(page('BLCA 最终候选区域 v3.4',body))

def report(s,cases):
    modified=s[s.changed.gt(0)]
    table='\n'.join(f'| {r.section_id} | {r.to_inside:,} | {r.to_outside:,} | {r.old_margin:,} | {r.new_margin:,} |' for r in modified.itertuples())
    no=', '.join(s[s.candidate_edges.eq(0)].section_id)
    holes=s[s.assigned_enclosed.gt(0)]
    text=f"""# 正式计算候选区域 v3.4

已完成63张切片、7个dataset、265,002个原始spot的最终交付。[浏览图形](index.html) · [完整逐点标签](final_spot_labels.tsv.gz) · [读取及合并示例](METHODS.md)。

相对 v3.3，**{len(modified)}张切片的{s.changed.sum():,}个点改变了候选位置**。其中{s.to_inside.sum():,}点移入 inside，{s.to_outside.sum():,}点移入 outside。47张非GSE319536切片和未通过配准的BC05、BC17均逐点核对。最终有{s.inside.sum():,}个inside、{s.outside.sum():,}个outside、{s.unlocalized.sum():,}个原有未定位点。

## 两处关键修订

BC18 复现了 array_row 14–36、array_col 23–57 的225点封闭蓝区。原作者图中心有116点落在Tumor Cell、19点落在Tumor Stroma。该区外邻域观测完整，所在作者肿瘤父区占多数，无正常粗区域反证。225点整体并入病灶内部区室，原组织成分仍单独保存。此前224点显示为200 μm外侧远区。最终225点均为候选内侧、距候选外缘超过200 μm，其中{cases['GSE319536_BC18']['exact']}点可给出精确负距离，余下仅能判定分带。[BC18局部图](by_dataset/GSE319536/GSE319536_BC18/revision.png)

BC13 的9、6、4、1、1点五个核心均保留。21点经真实组织邻接并入同一个{cases['GSE319536_BC13']['merged_n']:,}点ROI（{cases['GSE319536_BC13']['merged_roi']}），包括指定的两个单点barcode。原先引用这五个小ROI的106个200 μm外侧margin点，现在全部位于候选内侧且距新外缘超过200 μm。五个核心不再各自产生独立外缘。[BC13局部图](by_dataset/GSE319536/GSE319536_BC13/revision.png)

同一完整包围规则纳入的内部区室点数为 {', '.join(f'{r.section_id.split("_")[-1]} {r.assigned_enclosed}点' for r in holes.itertuples())}。BC14 有{cases['GSE319536_BC14']['single']}个 v3.3 新红点本身就是独立的单点观测组织块，全部保留inside，距离NA。B16 的 AGATAACACGGTATAA-1_4 保持outside。

## 修订切片与分带

| 切片 | 外→内 | 内→外 | v3.3 内外 margin 点数 | 最终内外 margin 点数 |
|---|---:|---:|---:|---:|
{table}

全队列可观测候选界面覆盖{s.candidate_edges.gt(0).sum()}张切片，{s.exact_distance.sum():,}点有精确候选距离，{s.band_200_determined.sum():,}点可判定200 μm分带。相对 v3.3，有{s.exited_margin.sum():,}点退出margin、{s.entered_margin.sum():,}点进入margin。无候选外缘切片为 {no}。在位置未改变的切片中，{s.loc[s.changed.eq(0),"withdrawn_outside_route"].sum():,}个外侧区室点因没有实际观测的外侧通路，撤回了旧外缘引用。缺测、截断和外侧通路不确定的局部区室继续保留位置，必要时距离NA。

## 使用边界

最终表适合作为全队列统一的计算候选区域。距离以当前候选ROI为条件，作者S1属于校准资料。作者粗区域边缘及非二元区室仍有定位不确定性，保留参考状态和原作者中心类别。T3 当前界面可能反映病灶内部组成／状态差异，尚不能等同于已确认的肿瘤—正常外缘。

本地交付包含逐片PNG/PDF、ROI和边段及7个dataset总览。GitHub写入接口返回403，当前Git环境没有可用写入凭据，本版交付为可上传目录。线上发布和线上链接尚未验证。
"""
    (OUT/'RESULTS.md').write_text(text)

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.targets=[]
    def handle_starttag(self,tag,attrs):
        self.targets.extend(v for k,v in attrs if k in ['href','src'])

def verify_links():
    targets=0
    for f in OUT.rglob('*.html'):
        parser=Links();parser.feed(f.read_text())
        for link in parser.targets:
            target=(f.parent/unquote(link)).resolve()
            assert target.is_relative_to(OUT) and target.is_file() and target.stat().st_size>0,(f,link)
            targets+=1
    print('Verified local HTML targets:',targets)

if __name__=='__main__':
    data=collect()
    stats=[]; cases={}
    for sec in M.index:
        s,c=validate(sec,data[data.section_id.eq(sec)])
        stats.append(s);cases[sec]=c
    summary=pd.DataFrame(stats)
    write_indexes(data,summary)
    report(summary,cases)
    verify_links()
    # Remove only transient shards produced by this final build. The master has
    # already been re-read and checked for original keys, coordinates and counts.
    keys=pd.read_csv(MASTER,sep='\t',usecols=['section_id','barcode'])
    assert keys.equals(data[['section_id','barcode']].reset_index(drop=True))
    for sec in M.index:
        (section_dir(sec)/'labels.tsv.gz').unlink(missing_ok=True)
    print(summary.select_dtypes('number').sum().to_dict())
    print('Cases:',{sec:c for sec,c in cases.items() if c})
