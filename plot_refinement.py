from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
ROOT=Path(__file__).resolve().parent;D=ROOT/'data';F=ROOT/'figures'
summary=pd.read_csv(D/'Refinement_Summary.csv').set_index('record')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':True,'axes.spines.right':True})
fig,axs=plt.subplots(2,2,figsize=(8.2,4.3))
a=axs[0,0]
a.plot([0,1],[.6110771885254435,summary.loc['BO06_finemesh','arithmetic_J_W']],'o-',c='red',label='BO06')
a.plot([0,1],[.6099833812668839,summary.loc['M201_finemesh','arithmetic_J_W']],'s-',c='blue',label='M201')
a.set_xticks([0,1],['Production','Fine']);a.set_ylabel('Arithmetic-mean J (W)');a.set_title('(a) Mesh refinement',loc='left');a.legend(frameon=False,fontsize=8)
a=axs[0,1];v=summary.loc[['BO06_finemesh','BO06_40wl'],'trapezoidal_J_W'].values
a.plot([20,40],v,'^-',c='blue');a.set_xticks([20,40]);a.set_xlabel('Wavelength samples');a.set_ylabel('Trapezoidal J (W)');a.set_title('(b) Same-mesh grid check: 0.4665%',loc='left')
for ax,col,label,panel in [(axs[1,0],'P_HgCdTe_W','HgCdTe power (W)','(c)'),(axs[1,1],'P_ITO_W','ITO power (W)','(d)')]:
 for record,color,legend in [('BO06_finemesh','red','20 wavelengths'),('BO06_40wl','blue','40 wavelengths')]:
  g=pd.read_csv(D/(record+'_merged.csv'));ax.plot(g.wavelength_um,g[col],c=color,lw=1.2,label=legend)
 ax.set_xlabel('Wavelength (µm)');ax.set_ylabel(label);ax.set_title(panel+' BO06 fine-mesh spectra',loc='left');ax.legend(frameon=False,fontsize=7,loc='lower left' if col=='P_HgCdTe_W' else 'upper right')
for ax in axs.flat:
 for sp in ax.spines.values():sp.set_visible(True)
 ax.tick_params(direction='out',which='both',top=True,right=True);ax.xaxis.set_minor_locator(AutoMinorLocator(2));ax.yaxis.set_minor_locator(AutoMinorLocator(2))
fig.tight_layout();fig.savefig(F/'refinement.png',dpi=300,bbox_inches='tight');fig.savefig(F/'refinement.pdf',bbox_inches='tight')
