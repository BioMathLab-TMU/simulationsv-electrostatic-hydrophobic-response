"""Recompute common-horizon structural and mobility tables from released inputs.

Run: python analysis/common84ms.py --output results/common84ms --verify
The original 28-ms release remains unchanged.
"""
from pathlib import Path
import argparse
import hashlib
import json
import h5py
import numpy as np
import pandas as pd
from scipy import stats
from interface_diagnostics import graph_metrics, self_test as graph_test
from mobility_diagnostic import fitted_mobility, self_test as mobility_test

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify',action='store_true')
    a=parser.parse_args()
    graph_test();mobility_test()
    manifest=json.loads((ROOT/'data/common84ms_availability.json').read_text())
    frames=[];run_rows=[];mob=[]
    for spec in manifest['runs']:
        path=ROOT/spec['path']
        with path.open('rb') as stream:
            assert hashlib.file_digest(stream,'sha256').hexdigest()==spec['sha256']
        identity={k:spec[k] for k in ('E','H','seed')}
        identity['condition']=f"E{spec['E']}H{spec['H']}"
        identity['run_id']=f"article_force_confirmatory_C30_E{spec['E']}_H{spec['H']}_seed_{spec['seed']}"
        local=[];xy=[]
        with h5py.File(path,'r') as f:
            names=sorted(f,key=lambda k:int(k[5:]))
            assert names==[f'step_{s}' for s in range(200,120001,200)]
            ids=f[names[0]]['ids'][()];types=f[names[0]]['types'][()]
            for name in names:
                g=f[name]
                np.testing.assert_array_equal(g['ids'][()],ids)
                np.testing.assert_array_equal(g['types'][()],types)
                p,r=g['positions'][()].astype(float),g['radii'][()].astype(float)
                assert np.isfinite(p).all()
                step=int(name[5:])
                for shell in (0,2,4,6):
                    row=dict(**identity,simulation_step=step,time_ms=step*.0007,shell_nm=shell,
                        **graph_metrics(p,r,types,shell=shell))
                    local.append(row)
                    if shell==2:frames.append(row)
                xy.append(p[types==1,:2]/1000)
        table=pd.DataFrame(local)
        metrics=list(graph_metrics(p,r,types))
        for shell,g in table.groupby('shell_nm'):
            g=g.sort_values('simulation_step')
            for horizon,count in ((28,200),(84,600)):
                tail=g.iloc[:count].iloc[-int(.3*count):]
                row=dict(**{k:identity[k] for k in ('E','H','condition','seed')},shell_nm=shell,horizon_ms=horizon,
                    late_frames=len(tail),late_start_ms=float(tail.time_ms.min()),late_stop_ms=float(tail.time_ms.max()))
                row.update({k:float(tail[k].mean()) for k in metrics})
                run_rows.append(row)
        xy=np.asarray(xy)
        for window,coords in (('full_28ms',xy[:200]),('late_28ms_60frames',xy[140:200]),
                              ('full_84ms',xy),('late_84ms_180frames',xy[-180:])):
            ds,lag,msd=fitted_mobility(coords,.00014)
            mob.append(dict(**identity,window=window,n_agents=len(ds),D_mean_um2_s=float(ds.mean()),
                D_median_um2_s=float(np.median(ds)),negative_slope_fraction=float(np.mean(ds<0))))
        print('RECOMPUTED',identity['condition'],spec['seed'],flush=True)
    runs=pd.DataFrame(run_rows)
    primary=runs.loc[runs.shell_nm.eq(2)]
    contrasts=[]
    for horizon,g in primary.groupby('horizon_ms'):
        for metric in metrics:
            pv=g.pivot(index='seed',columns=['E','H'],values=metric)
            effects={'H_at_E0':pv[0,1]-pv[0,0],'H_at_E1':pv[1,1]-pv[1,0],
                'E_at_H0':pv[1,0]-pv[0,0],'E_at_H1':pv[1,1]-pv[0,1],
                'E_H_interaction':pv[1,1]-pv[1,0]-pv[0,1]+pv[0,0]}
            for name,v in effects.items():
                row=dict(horizon_ms=horizon,metric=metric,contrast=name,mean=float(v.mean()),sd=float(v.std(ddof=1)))
                row.update({f'seed_{s}':float(v.loc[s]) for s in (9310002,9310003,9310004)})
                width=stats.t.ppf(.975,2)*v.std(ddof=1)/np.sqrt(3)
                row.update(ci95_low=float(v.mean()-width),ci95_high=float(v.mean()+width))
                contrasts.append(row)
    outputs={'frame_metrics.csv':pd.DataFrame(frames),'run_metrics.csv':runs,
        'run_level_mobility.csv':pd.DataFrame(mob),'paired_contrasts.csv':pd.DataFrame(contrasts)}
    a.output.mkdir(parents=True,exist_ok=True)
    for name,table in outputs.items():
        if a.verify:
            old=pd.read_csv(ROOT/'data/processed/common84ms'/name)
            pd.testing.assert_frame_equal(table.reset_index(drop=True)[old.columns],old,
                check_dtype=False,atol=1e-12,rtol=1e-12)
        table.to_csv(a.output/name,index=False)
    print('COMMON84MS_REPRODUCTION_PASSED verify='+str(a.verify))


if __name__=='__main__':main()
