from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
ROOT=Path(__file__).resolve().parent
c=pd.read_csv(ROOT/'data/Corrected_Candidate_Spectra_480.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelweight':'bold','axes.linewidth':1,'axes.spines.top':True,'axes.spines.right':True})
fig,axs=plt.subplots(1,2,figsize=(8.2,3.65));fig.subplots_adjust(left=.075,right=.985,bottom=.16,top=.97,wspace=.30)
for ax in axs:
 ax.tick_params(which='both',direction='out',top=True,right=True,labelsize=8)
 ax.tick_params(which='major',length=4,width=1);ax.tick_params(which='minor',length=2.5,width=.8)
 ax.xaxis.set_major_locator(MultipleLocator(.5));ax.xaxis.set_minor_locator(AutoMinorLocator(5));ax.yaxis.set_minor_locator(AutoMinorLocator(5))
 ax.set_xlim(2.47,5.03);ax.set_xlabel('Wavelength (µm)',fontsize=8.5)
 for s in ax.spines.values():s.set_visible(True)
def inset(ax,labels,vals,pos):
 a=ax.inset_axes(pos);a.plot(range(len(vals)),vals,'^-',color='blue',lw=1,ms=3);a.set_xticks(range(len(vals)),labels,fontsize=6)
 a.tick_params(which='both',direction='inout',top=True,right=True,labelsize=6,length=3);a.set_ylabel('AUC (W µm)',fontsize=6,labelpad=2)
 a.margins(x=.1,y=.16)
 for s in a.spines.values():s.set_visible(True)
vals=[]
for key,color,label in [('BO06','red','BO06'),('RS03','blue','RS03')]:
 g=c[c.candidate_id==key].sort_values('wavelength_um');axs[0].plot(g.wavelength_um,g.P_HgCdTe_W,c=color,lw=1.3,label=label);vals.append(np.trapezoid(g.P_HgCdTe_W,g.wavelength_um))
axs[0].set_ylim(.21,.80);axs[0].set_ylabel('HgCdTe absorbed power (W)',fontsize=8.5)
inset(axs[0],['BO06','RS03'],vals,(.13,.24,.33,.26));axs[0].legend(loc='lower right',bbox_to_anchor=(.87,.16),frameon=False,fontsize=7,handlelength=1.8)
g=c[c.candidate_id=='BO06'].sort_values('wavelength_um');vals=[]
for col,color,label in [('P_HgCdTe_W','red','HgCdTe'),('P_ITO_W','blue','ITO'),('A_total','black','Total')]:
 axs[1].plot(g.wavelength_um,g[col],color=color,lw=1.3,label=label);vals.append(np.trapezoid(g[col],g.wavelength_um))
axs[1].set_ylim(-.01,.88);axs[1].set_ylabel('Absorbed power (W)',fontsize=8.5)
inset(axs[1],['HgCdTe','ITO','Total'],vals,(.13,.25,.36,.26));axs[1].legend(loc='lower right',bbox_to_anchor=(.85,.24),frameon=False,fontsize=7,handlelength=1.8)
fig.savefig(ROOT/'figures/paired_spectra.png',dpi=320,bbox_inches='tight');fig.savefig(ROOT/'figures/paired_spectra.pdf',bbox_inches='tight')
