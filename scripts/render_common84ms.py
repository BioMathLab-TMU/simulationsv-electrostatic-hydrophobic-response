"""Reproduce article-1 figures at the common 84-ms observation horizon."""
from pathlib import Path
import argparse
import json
import os
import tempfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
import numpy as np
import pandas as pd

CONDITIONS=[(0,0),(0,1),(1,0),(1,1)]
COLORS=['#737373','#E76F51','#2A9D8F','#195E56']
MARKERS=['o','s','^','D']
STYLES=[':','--','-.','-']


def finish(fig,out,name):
    fig.canvas.draw()
    regular=Path('C:/Windows/Fonts/arial.ttf')
    if regular.exists():
        for item in fig.findobj(match=Text):
            file=regular.with_name('arialbd.ttf') if item.get_fontweight()=='bold' else regular
            item.set_fontproperties(font_manager.FontProperties(fname=str(file),size=item.get_fontsize()))
    out.mkdir(parents=True,exist_ok=True)
    for ext in ('pdf','png','tif'):
        fd,temp=tempfile.mkstemp(prefix='.figure84_',suffix='.'+ext,dir=out)
        os.close(fd)
        try:
            options={'dpi':600 if ext=='tif' else 180}
            if ext=='tif': options['pil_kwargs']={'compression':'tiff_lzw'}
            fig.savefig(temp,**options)
            os.replace(temp,out/f'{name}.{ext}')
        finally:
            if os.path.exists(temp): os.unlink(temp)
    plt.close(fig)


def style(ax,panel):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.15)
    ax.set_title(panel,loc='left',fontweight='bold',fontsize=12)


def dots(ax,table,column):
    for x,((e,h),color,marker) in enumerate(zip(CONDITIONS,COLORS,MARKERS)):
        values=table.loc[table.E.eq(e)&table.H.eq(h),column].to_numpy(float)
        assert len(values)==3
        ax.scatter(x+np.linspace(-.08,.08,3),values,s=30,facecolor='white',edgecolor=color,marker=marker,zorder=4)
        ax.errorbar(x,values.mean(),yerr=values.std(ddof=1),color=color,marker='_',markersize=17,capsize=3,lw=1.3,zorder=3)
    ax.set_xticks(range(4),[f'E{e}H{h}' for e,h in CONDITIONS])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,
        'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42})
    frames=pd.read_csv(a.results/'frame_metrics.csv')
    runs=pd.read_csv(a.results/'run_metrics.csv')
    runs=runs.loc[runs.horizon_ms.eq(84)&runs.shell_nm.eq(2)]
    assert len(runs)==12 and len(frames)==7200
    metrics=['largest_component_fraction','syn_nonisolated_fraction','largest_component_sv_fraction']
    labels=['Largest component /\nall agents','Non-isolated Syn /\nall Syn','SVs in largest component /\nall SVs']
    fig,axes=plt.subplots(1,3,figsize=(7.5,3.25),layout='constrained')
    for ax,key,label,panel in zip(axes,metrics,labels,'ABC'):
        dots(ax,runs,key); ax.set(ylabel=label,ylim=(0,1)); style(ax,panel)
    finish(fig,a.output,'force_factorial_endpoints')

    fig=plt.figure(figsize=(7.5,7.0),layout='constrained')
    grid=fig.add_gridspec(2,2)
    axa=fig.add_subplot(grid[0,0]); axb=fig.add_subplot(grid[0,1]); axc=fig.add_subplot(grid[1,:])
    series={}
    for (e,h),color,marker,ls in zip(CONDITIONS,COLORS,MARKERS,STYLES):
        g=frames.loc[frames.E.eq(e)&frames.H.eq(h)]
        s=g.groupby('simulation_step').agg(time=('time_ms','mean'),mean=('largest_component_fraction','mean'),sd=('largest_component_fraction','std'))
        series[e,h]=s
        axa.plot(s.time,s['mean'],color=color,ls=ls,lw=1.4,label=f'E{e}H{h}')
        axa.fill_between(s.time,s['mean']-s.sd,s['mean']+s.sd,color=color,alpha=.12)
        v=runs.loc[runs.E.eq(e)&runs.H.eq(h)]
        x=v.syn_nonisolated_fraction;y=v.largest_component_sv_fraction
        axb.scatter(x,y,marker=marker,edgecolor=color,facecolor='white',s=31,zorder=4)
        axb.errorbar(x.mean(),y.mean(),xerr=x.std(ddof=1),yerr=y.std(ddof=1),marker=marker,color=color,markersize=5,capsize=3)
    ref=series[1,0]['mean']+series[0,1]['mean']-series[0,0]['mean']
    axa.plot(series[0,0].time,ref,color='black',ls=(0,(4,2)),lw=1,label='Additive reference')
    axa.axvspan(58.94,84,color='gray',alpha=.08,zorder=0)
    axa.set(xlabel='Time (ms)',ylabel='Largest component / all agents',xlim=(0,84),ylim=(0,1),xticks=[0,28,56,84])
    axa.legend(loc='upper left',frameon=False)
    axb.set(xlabel='Non-isolated Syn / all Syn',ylabel='SVs in largest component / all SVs',xlim=(.25,1),ylim=(0,1))
    joint=frames.loc[frames.E.eq(1)&frames.H.eq(1)]
    for (seed,g),ls in zip(joint.groupby('seed'),('-','--',':')):
        axc.plot(g.time_ms,g.largest_component_fraction,color='#737373',ls=ls,lw=.85,alpha=.7,label=f'Seed {seed}')
    s=series[1,1]
    axc.plot(s.time,s['mean'],color=COLORS[3],lw=1.65,label='Mean')
    axc.axvspan(58.94,84,color='gray',alpha=.08,zorder=0)
    axc.axvline(28,color='gray',ls='-.',lw=.85)
    axc.set(xlabel='Time (ms)',ylabel='Largest component / all agents',xlim=(0,84),ylim=(0,1),xticks=[0,14,28,42,56,70,84])
    axc.legend(loc='upper left',frameon=False,ncol=2)
    for ax,panel in zip((axa,axb,axc),'ABC'):style(ax,panel)
    finish(fig,a.output,'force_factorial_timecourse')

    fig,axes=plt.subplots(1,3,figsize=(7.5,3.25),layout='constrained')
    columns=['syn_near_sv_fraction','syn_multivesicle_neighbor_fraction','edges_SV']
    labels=['Syn with at least one SV\nneighbor / all Syn','Syn with at least two SV\nneighbors / all Syn','Syn-SV edges\nper snapshot']
    for k,(ax,column,label) in enumerate(zip(axes,columns,labels)):
        dots(ax,runs,column)
        ax.set(ylabel=label,ylim=(0,1) if k<2 else (0,float(runs.edges_SV.max())*1.15))
        style(ax,'ABC'[k])
    finish(fig,a.output,'force_interface_diagnostics')

    mob=pd.read_csv(a.results/'run_level_mobility.csv')
    full=mob.loc[mob.window.eq('full_84ms')]
    fig,ax=plt.subplots(figsize=(6.6,3.6),layout='constrained')
    dots(ax,full,'D_mean_um2_s')
    reference=json.loads(a.reference.read_text())
    ax.axhline(reference['with_SV']['mean'],color='black',ls='--',lw=1.1,label='Experiment with SVs (mean)')
    nominal=1.380649e-23*310.15/(6*np.pi*.35*3.972e-9)*1e12
    ax.axhline(nominal,color='gray',ls=':',lw=1.1,label='Nominal free Syn input')
    ax.set(ylabel='Apparent short-lag mobility ($\\mu$m$^2$/s)',ylim=(0,.185))
    ax.legend(loc='upper right',frameon=False);style(ax,'')
    finish(fig,a.output,'force_mobility_external_comparison')
    provenance=dict(comparison_horizon_ms=84,independent_seeds=3,late_window_ms=[58.94,84],
        late_frames=180,errors='between-seed SD',tiff_dpi=600,compression='LZW',
        palette={f'E{e}H{h}':c for (e,h),c in zip(CONDITIONS,COLORS)},smoothing='none')
    temp=a.output/'figure_provenance.json.tmp'
    temp.write_text(json.dumps(provenance,indent=2)+'\n')
    os.replace(temp,a.output/'figure_provenance.json')
    print('FOUR_COMMON84MS_FIGURES_COMPLETE')


if __name__=='__main__':main()
