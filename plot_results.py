from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
ROOT=Path(__file__).resolve().parent;D=ROOT/'data';F=ROOT/'figures'
summary=json.loads((ROOT/'analysis_summary.json').read_text())
bench=pd.read_csv(D/'Final_Model_Benchmarks.csv');lc=pd.read_csv(D/'Learning_Curve_Runs.csv')
bp=pd.read_csv(D/'Final_Blind_Predictions.csv');imp=pd.read_csv(D/'Validation_Permutation_Importance.csv')
cand=pd.read_csv(D/'Candidate_Spectra_With_Audit_Predictions.csv');res=pd.read_csv(D/'Corrected_All_24_Candidate_Results.csv')
designs=pd.read_csv(D/'Main300_Final_Geometry_Partitions.csv')
st=summary['statistics']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':240})
teal='#12686e';gold='#be7c29';blue='#3874a3';gray='#7f8d91'
names={'P_HgCdTe_W':'HgCdTe','P_ITO_W':'ITO'}
labels={'ito_w_um':'ITO lateral size','ito_h_um':'ITO height','cell_w_um':'Cell width','cell_d_um':'Cell depth','wavelength_um':'Wavelength'}
def save(fig,name):
 fig.tight_layout();fig.savefig(F/(name+'.png'),bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(7.7,3.4))
for target,color in zip(names,[teal,gold]):
 g=lc[lc.target==target].groupby('training_designs').MAE.agg(['mean','std'])
 ax.errorbar(g.index[:-1],g['mean'].iloc[:-1],yerr=g['std'].iloc[:-1],fmt='o-',color=color,capsize=4,label=names[target])
 ax.plot(g.index[-2:],g['mean'].iloc[-2:],'o-',color=color)
ax.set_xticks([30,75,125,200,210]);ax.set_xlabel('Independent training geometries');ax.set_ylabel('Validation MAE (W)');ax.legend(frameon=False);ax.grid(alpha=.18)
save(fig,'learning')
fig,axs=plt.subplots(1,2,figsize=(8.2,3.8))
for ax,target,color in zip(axs,names,[teal,gold]):
 g=bp[bp.target==target];ax.scatter(g.actual,g.predicted,s=7,alpha=.5,c=color,edgecolors='none')
 low=min(g.actual.min(),g.predicted.min());high=max(g.actual.max(),g.predicted.max())
 ax.plot([low,high],[low,high],'--',c=gray,lw=1)
 row=bench[(bench.target==target)&(bench.model==summary['selected_models'][target])&(bench.split=='blind')].iloc[0]
 ax.text(.04,.96,f"MAE {row.MAE:.5f} W\nRMSE {row.RMSE:.5f} W\nR² {row.R2:.5f}",transform=ax.transAxes,va='top',fontsize=8.5)
 ax.set_title(names[target]);ax.set_xlabel('COMSOL power (W)');ax.set_ylabel('Predicted power (W)');ax.set_aspect('equal',adjustable='box')
save(fig,'blind')
fig,axs=plt.subplots(1,2,figsize=(8.3,3.5))
ax=axs[0];x=np.arange(2)
for off,model,color in [(-.18,'Wavelength-only Extra Trees',gray),(.18,'Extra Trees',teal)]:
 vals=[bench[(bench.target==t)&(bench.model==(model if off<0 else summary['selected_models'][t]))&(bench.split=='blind')].MAE.iloc[0] for t in names]
 ax.bar(x+off,vals,.34,label='Wavelength only' if off<0 else 'Full features',color=color)
ax.set_xticks(x,list(names.values()));ax.set_ylabel('Blind MAE (W)');ax.set_title('(a) Incremental geometry information');ax.legend(fontsize=8,frameon=False)
g=imp[imp.target=='P_HgCdTe_W'].sort_values('MAE_increase_mean')
axs[1].barh([labels[f] for f in g.feature],g.MAE_increase_mean,xerr=g.MAE_increase_sd,color=teal,capsize=3)
axs[1].set_xlabel('Validation MAE increase (W)');axs[1].set_title('(b) HgCdTe permutation importance')
save(fig,'importance')
fig,axs=plt.subplots(1,2,figsize=(8.2,3.65),gridspec_kw={'width_ratios':[1.2,1]})
ax=axs[0]
for pos,prefix,color in [(0,'BO',teal),(1,'RS',gold)]:
 g=res[res.candidate_id.str.startswith(prefix)].sort_values('candidate_id')
 jitter=np.linspace(-.10,.10,len(g))
 ax.scatter(pos+jitter,g.verified_J,s=24,c=color,zorder=3)
 ax.hlines(g.verified_J.mean(),pos-.23,pos+.23,colors='black',lw=1.2)
ax.set_xticks([0,1],['Guided (n = 12)','Random (n = 12)']);ax.set_ylabel('Verified arithmetic-mean J (W)');ax.grid(axis='y',alpha=.2);ax.set_title('(a) Individual verification results')
lo,hi=st['CI'];mean=st['difference']
axs[1].errorbar([mean],[0],xerr=[[mean-lo],[hi-mean]],fmt='o',color=teal,capsize=5)
axs[1].axvline(0,ls='--',c=gray,lw=1);axs[1].set_yticks([]);axs[1].set_xlim(-.001,.023);axs[1].set_ylim(-1,1)
axs[1].text(mean,.25,f"{mean:.6f} W",ha='center',fontsize=9)
axs[1].text(mean,-.33,f"95% CI {lo:.6f}–{hi:.6f}",ha='center',fontsize=7.6)
axs[1].set_xlabel('Guided minus random J (W)');axs[1].set_title('(b) Bootstrap mean difference')
save(fig,'distribution')
fig,ax=plt.subplots(figsize=(7.7,3.7))
for key,color,ls in [('BO06',teal,'-'),(st['best_random'],gold,'--')]:
 g=cand[cand.candidate_id==key].sort_values('wavelength_um')
 ax.plot(g.wavelength_um,g.P_HgCdTe_W,ls,color=color,label=key,lw=1.6)
ax.set_xlabel('Wavelength (µm)');ax.set_ylabel('HgCdTe power (W)');ax.legend(frameon=False);ax.grid(alpha=.15)
save(fig,'best_spectra')
fig,ax=plt.subplots(figsize=(7.7,3.7))
g=cand[cand.candidate_id=='BO06'].sort_values('wavelength_um')
for col,color,ls,label in [('P_HgCdTe_W',teal,'-','HgCdTe'),('P_ITO_W',gold,'-','ITO'),('A_total',gray,'--','Total absorptance × 1 W')]:
 ax.plot(g.wavelength_um,g[col],ls,color=color,label=label,lw=1.5)
ax.set_xlabel('Wavelength (µm)');ax.set_ylabel('Absorbed power / equivalent power (W)');ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.15)
save(fig,'material_spectra')
# Supplementary coverage of the independent four-dimensional geometry.
cols=['ito_w_um','ito_h_um','cell_w_um','cell_d_um']
fig,axs=plt.subplots(4,4,figsize=(8.3,8.3))
for i,y in enumerate(cols):
 for j,x in enumerate(cols):
  ax=axs[i,j]
  if i==j:ax.hist(designs[x],bins=14,color=teal,alpha=.8)
  else:ax.scatter(designs[x],designs[y],s=4,alpha=.45,c=teal)
  if i==3:ax.set_xlabel(labels[x],fontsize=8)
  if j==0:ax.set_ylabel(labels[y],fontsize=8)
  ax.tick_params(labelsize=7)
save(fig,'s_coverage')
fig,axs=plt.subplots(2,2,figsize=(8.1,5.6))
for ax,col in zip(axs.flat,cols):
 for split,color in [('train',teal),('validation',gold),('blind_test',blue)]:
  vals=np.sort(designs.loc[designs.split==split,col]);ax.step(vals,np.arange(1,len(vals)+1)/len(vals),where='post',c=color,label=split)
 ax.set_xlabel(labels[col]+' (µm)');ax.set_ylabel('Empirical CDF');ax.grid(alpha=.15)
axs[0,0].legend(fontsize=7,frameon=False)
save(fig,'s_partitions')
fig,axs=plt.subplots(1,2,figsize=(8.3,4))
order=bench.model.drop_duplicates().tolist()
for ax,t in zip(axs,names):
 g=bench[bench.target==t]
 for off,split,col in [(-.17,'validation',teal),(.17,'blind',gold)]:
  vals=[g[(g.model==m)&(g.split==split)].MAE.iloc[0] for m in order]
  ax.barh(np.arange(len(order))+off,vals,.32,color=col,label=split)
 ax.set_yticks(np.arange(len(order)),['Mean','Polynomial Ridge','Random Forest','Extra Trees','Hist. boosting','Wavelength only'],fontsize=8)
 ax.set_xlabel('MAE (W)');ax.set_title(names[t]);ax.invert_yaxis();ax.legend(frameon=False,fontsize=8)
save(fig,'s_benchmarks')
fig,axs=plt.subplots(1,2,figsize=(8.2,3.6))
for ax,t,col in zip(axs,names,[teal,gold]):
 g=bp[bp.target==t];ax.scatter(g.wavelength_um,g.signed_error,s=5,alpha=.35,c=col,edgecolors='none')
 ax.axhline(0,c=gray,ls='--',lw=1);ax.set_xlabel('Wavelength (µm)');ax.set_ylabel('Predicted − COMSOL (W)');ax.set_title(names[t])
save(fig,'s_residuals')
for prefix,filename in [('BO','s_guided'),('RS','s_random')]:
 fig,axs=plt.subplots(1,2,figsize=(8.3,3.7))
 for key,g in cand[cand.candidate_id.str.startswith(prefix)].groupby('candidate_id'):
  g=g.sort_values('wavelength_um');axs[0].plot(g.wavelength_um,g.P_HgCdTe_W,label=key,lw=.9);axs[1].plot(g.wavelength_um,g.P_ITO_W,label=key,lw=.9)
 for ax,label in zip(axs,['HgCdTe power (W)','ITO power (W)']):
  ax.set_xlabel('Wavelength (µm)');ax.set_ylabel(label);ax.grid(alpha=.15)
 axs[1].legend(ncol=3,fontsize=6,frameon=False,loc='upper center',bbox_to_anchor=(.45,-.20))
 save(fig,filename)
fig,axs=plt.subplots(1,2,figsize=(8.2,3.5))
for ax,t,c in zip(axs,names,[teal,gold]):
 g=lc[lc.target==t]
 ax.scatter(g.training_designs,g.MAE,color=c,label='Individual subset fits')
 means=g.groupby('training_designs').MAE.mean();ax.plot(means.index,means.values,color=c)
 ax.set_xlabel('Independent training geometries');ax.set_ylabel('Validation MAE (W)');ax.set_title(names[t]);ax.set_xticks([30,75,125,200,210]);ax.grid(alpha=.15)
save(fig,'s_learning_replicates')
fig,axs=plt.subplots(1,2,figsize=(8.2,3.7))
for ax,col,actual,label in [(axs[0],'pred_mean_HgCdTe_W','mean_HgCdTe','Mean HgCdTe power'),(axs[1],'pred_objective_J','verified_J','Mean objective J')]:
 for prefix,c in [('BO',teal),('RS',gold)]:
  g=res[res.candidate_id.str.startswith(prefix)];ax.scatter(g[actual],g[col],label=prefix,c=c,s=22)
 a=min(res[actual].min(),res[col].min());b=max(res[actual].max(),res[col].max());ax.plot([a,b],[a,b],'--',c=gray,lw=1)
 ax.set_xlabel('Verified (W)');ax.set_ylabel('Original frozen prediction (W)');ax.set_title(label);ax.legend(frameon=False,fontsize=8)
save(fig,'s_candidate_errors')
print('Generated figures',len(list(F.glob('*.png'))))
