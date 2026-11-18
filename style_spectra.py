"""Reference-matched spectral style; numerical curves remain unchanged."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator,MultipleLocator
ROOT=Path(__file__).resolve().parent
D=ROOT/'data';F=ROOT/'figures'
c=pd.read_csv(D/'Corrected_Candidate_Spectra_480.csv')
r=pd.read_csv(D/'Corrected_All_24_Candidate_Results.csv')
best=r[r.candidate_id.str.startswith('RS')].sort_values('verified_J').iloc[-1].candidate_id
palette=['#ff0000','#ff6600','#ffa500','#e6d800','#65e600','#00b4e6','#0033dd','#9b00c9','#db2386','#7a471c','#777777','#000000']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelweight':'bold','axes.linewidth':1.4,'savefig.dpi':300,'axes.spines.top':True,'axes.spines.right':True})
areas=[]
def style(ax,ylabel):
 for s in ax.spines.values():s.set_visible(True);s.set_color('black');s.set_linewidth(1.4)
 ax.tick_params(which='both',direction='out',top=True,right=True,labeltop=False,labelright=False)
 ax.tick_params(which='major',length=7,width=1.3,pad=5)
 ax.tick_params(which='minor',length=4,width=1.0)
 ax.xaxis.set_major_locator(MultipleLocator(.5));ax.xaxis.set_minor_locator(AutoMinorLocator(5));ax.yaxis.set_minor_locator(AutoMinorLocator(5))
 ax.set_xlim(2.47,5.03);ax.set_xlabel('Wavelength (µm)');ax.set_ylabel(ylabel);ax.grid(False)
def auc(g,col):
 return float(np.trapezoid(g[col],g.wavelength_um))
def inset(ax,labels,values,position=(.12,.16,.29,.28)):
 ins=ax.inset_axes(position)
 ins.plot(range(len(values)),values,color='#0033ee',marker='^',lw=1.5,ms=4)
 ins.set_xticks(range(len(values)),labels,rotation=45 if len(values)>3 else 0)
 ins.set_ylabel('AUC (W µm)',fontsize=8,labelpad=2)
 ins.tick_params(which='both',direction='inout',top=True,right=True,labelsize=7,length=5,width=1)
 ins.yaxis.set_minor_locator(AutoMinorLocator(2))
 for s in ins.spines.values():s.set_visible(True);s.set_linewidth(1.1)
 ins.margins(x=.05,y=.18)
 ins.ticklabel_format(axis='y',style='plain',useOffset=False)
 return ins
def save(fig,key):
 fig.savefig(F/(key+'.png'),bbox_inches='tight');fig.savefig(F/(key+'.pdf'),bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(7.6,5.25));fig.subplots_adjust(left=.13,right=.97,bottom=.14,top=.97)
vals=[]
for key,color,label in [('BO06','#ff0000','BO06 (best guided)'),(best,'#0033dd',best+' (best random)')]:
 g=c[c.candidate_id==key].sort_values('wavelength_um');ax.plot(g.wavelength_um,g.P_HgCdTe_W,c=color,lw=1.8,label=label);vals.append(auc(g,'P_HgCdTe_W'))
style(ax,'HgCdTe absorbed power (W)');ax.set_ylim(.21,.80)
inset(ax,['BO06',best],vals)
ax.legend(loc='lower right',bbox_to_anchor=(.93,.12),frameon=False,fontsize=10)
save(fig,'best_spectra')
fig,ax=plt.subplots(figsize=(7.6,5.25));fig.subplots_adjust(left=.13,right=.97,bottom=.14,top=.97)
g=c[c.candidate_id=='BO06'].sort_values('wavelength_um');vals=[]
for col,color,label in [('P_HgCdTe_W','#ff0000','HgCdTe'),('P_ITO_W','#0033dd','ITO'),('A_total','#000000','Total absorbed power')]:
 ax.plot(g.wavelength_um,g[col],c=color,lw=1.8,label=label);vals.append(auc(g,col))
style(ax,'Absorbed power (W)');ax.set_ylim(-.015,.88)
inset(ax,['HgCdTe','ITO','Total'],vals,position=(.12,.24,.29,.28))
ax.legend(loc='lower right',bbox_to_anchor=(.94,.25),frameon=False,fontsize=10)
save(fig,'material_spectra')
for prefix,key in [('BO','s_guided'),('RS','s_random')]:
 fig,axs=plt.subplots(2,1,figsize=(7.6,9.2));fig.subplots_adjust(left=.14,right=.97,bottom=.07,top=.97,hspace=.30)
 for ax,col,ylab,panel in zip(axs,['P_HgCdTe_W','P_ITO_W'],['HgCdTe absorbed power (W)','ITO absorbed power (W)'],['(a)','(b)']):
  ids=[];values=[]
  for color,(cid,g) in zip(palette,c[c.candidate_id.str.startswith(prefix)].groupby('candidate_id',sort=True)):
   g=g.sort_values('wavelength_um');ax.plot(g.wavelength_um,g[col],c=color,lw=1.35,label=cid);ids.append(cid[-2:]);values.append(auc(g,col));areas.append({'candidate_id':cid,'material':col,'trapezoidal_spectral_AUC_W_um':values[-1]})
  style(ax,ylab);ax.text(.02,.96,panel,transform=ax.transAxes,va='top',fontweight='bold')
  if col=='P_HgCdTe_W':
   ax.set_ylim(.20,.81);pos=(.11,.18,.30,.27);anchor=(.82,.14);loc='lower right'
  else:
   ax.set_ylim(-.003,max(c.loc[c.candidate_id.str.startswith(prefix),col])*1.15);pos=(.35,.68,.30,.23);anchor=(.96,.91);loc='upper right'
  ii=inset(ax,ids,values,position=pos);ii.set_xlabel('Candidate number',fontsize=8,labelpad=1)
  ii.set_xticks([0,3,6,9,11],[ids[j] for j in [0,3,6,9,11]],rotation=0)
  ax.legend(loc=loc,bbox_to_anchor=anchor,ncol=2,frameon=False,fontsize=8.2,handlelength=2.1,columnspacing=.8,labelspacing=.35)
 save(fig,key)
pd.DataFrame(areas).to_csv(D/'Candidate_Spectral_AUC.csv',index=False)
print('Restyled 4 spectral figures with boxed axes, reference palette, internal legends and computed AUC insets.')
