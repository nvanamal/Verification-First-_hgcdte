"""Compact horizontal supplementary spectra, retaining all 24 candidates."""
from pathlib import Path
p=Path(__file__).resolve().parent
source=(p/'style_spectra.py').read_text()
# Reuse source-data loading, four-sided axis styling, palette and AUC helpers.
exec(source[:source.index('fig,ax=plt.subplots')])
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7})
for prefix,key in [('BO','s_guided_compact'),('RS','s_random_compact')]:
 fig,axs=plt.subplots(1,2,figsize=(8.2,3.7));fig.subplots_adjust(left=.08,right=.985,bottom=.17,top=.96,wspace=.32)
 for ax,col,label,panel in zip(axs,['P_HgCdTe_W','P_ITO_W'],['HgCdTe absorbed power (W)','ITO absorbed power (W)'],['(a)','(b)']):
  ids=[];values=[]
  for color,(cid,g) in zip(palette,c[c.candidate_id.str.startswith(prefix)].groupby('candidate_id',sort=True)):
   g=g.sort_values('wavelength_um');ax.plot(g.wavelength_um,g[col],c=color,lw=1,label=cid);ids.append(cid[-2:]);values.append(auc(g,col))
  style(ax,label);ax.tick_params(labelsize=7,length=4);ax.tick_params(which='minor',length=2);ax.text(.03,.96,panel,transform=ax.transAxes,va='top',fontsize=8)
  if col=='P_HgCdTe_W':
   ax.set_ylim(.20,.81);pos=(.12,.23,.29,.25);anchor=(.85,.14);loc='lower right'
  else:
   ax.set_ylim(-.003,max(c.loc[c.candidate_id.str.startswith(prefix),col])*1.15);pos=(.22,.67,.30,.22);anchor=(.98,.92);loc='upper right'
  ii=inset(ax,ids,values,position=pos);ii.set_xlabel('Candidate',fontsize=5.5,labelpad=0);ii.set_ylabel('AUC (W µm)',fontsize=5.5,labelpad=0)
  ii.tick_params(labelsize=5.5,length=3);ii.set_xticks([0,5,11],['01','06','12'],rotation=0)
  ax.legend(loc=loc,bbox_to_anchor=anchor,ncol=2,frameon=False,fontsize=5.8,handlelength=1.4,columnspacing=.5,labelspacing=.25)
 save(fig,key)
