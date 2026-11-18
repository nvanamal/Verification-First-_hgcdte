from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import spearmanr,kendalltau,mannwhitneyu
ROOT=Path(__file__).resolve().parent;D=ROOT/'data';seed=2026
c=pd.read_csv(D/'Corrected_All_24_Candidate_Results.csv');m=pd.read_csv(D/'Main300_Design_Mean_Powers.csv').set_index('design_id').loc['M201'];b=pd.read_csv(D/'Blind_Geometry_Objectives.csv');spectra=pd.read_csv(D/'Corrected_Candidate_Spectra_480.csv')
rows=[]
for alpha in [0,.5,1,2]:
 j=c.mean_HgCdTe-alpha*c.mean_ITO;g=j[c.candidate_id.str.startswith('BO')];r=j[c.candidate_id.str.startswith('RS')];v=c.loc[j.idxmax(),'candidate_id'];bo=c[c.candidate_id=='BO06'].iloc[0]
 rows.append({'alpha':alpha,'best_among_24':v,'guided_mean_J':g.mean(),'random_mean_J':r.mean(),'guided_minus_random':g.mean()-r.mean(),'all_144_pairs_favor_guided':bool(g.min()>r.max()),'BO06_J':bo.mean_HgCdTe-alpha*bo.mean_ITO,'M201_J':m.P_HgCdTe_W-alpha*m.P_ITO_W,'BO06_minus_M201':bo.mean_HgCdTe-alpha*bo.mean_ITO-(m.P_HgCdTe_W-alpha*m.P_ITO_W)})
pd.DataFrame(rows).to_csv(D/'Objective_Weight_Sensitivity.csv',index=False)
trap=[]
for k,g in spectra.groupby('candidate_id'):
 g=g.sort_values('wavelength_um');trap.append({'candidate_id':k,'trapezoidal_J':np.trapezoid(g.P_HgCdTe_W-g.P_ITO_W,g.wavelength_um)/2.5})
t=c[['candidate_id','verified_J']].merge(pd.DataFrame(trap),on='candidate_id');t['arithmetic_rank']=t.verified_J.rank(ascending=False).astype(int);t['trapezoidal_rank']=t.trapezoidal_J.rank(ascending=False).astype(int);t.to_csv(D/'Arithmetic_Trapezoidal_Rankings.csv',index=False)
a=c[c.candidate_id.str.startswith('BO')].verified_J.to_numpy();r=c[c.candidate_id.str.startswith('RS')].verified_J.to_numpy()
main=pd.read_csv(D/'Main300_Final_Geometry_Partitions.csv');geo=['ito_w_um','ito_h_um','cell_w_um','cell_d_um'];dup=0
for row in c[geo].to_numpy():dup+=int(np.isclose(main[geo].to_numpy(),row,rtol=0,atol=1e-12).all(axis=1).any())
k=int(np.ceil(.1*len(b)));actual=set(b.nlargest(k,'actual_J').design_id);pred=set(b.nlargest(k,'predicted_J').design_id)
rng=np.random.default_rng(seed);errors=abs(b.J_error.to_numpy());idx=rng.integers(0,len(errors),(100000,len(errors)));maes=errors[idx].mean(axis=1)
stats={'blind_geometries':len(b),'blind_J_MAE':float(errors.mean()),'blind_J_RMSE':float(np.sqrt(np.mean(b.J_error**2))),'blind_J_spearman':float(spearmanr(b.actual_J,b.predicted_J).statistic),'blind_J_kendall':float(kendalltau(b.actual_J,b.predicted_J).statistic),'top_decile_k':k,'top_decile_overlap':len(actual&pred),'top_decile_precision':len(actual&pred)/k,'blind_J_MAE_bootstrap_95CI':np.quantile(maes,[.025,.975]).tolist(),'candidate_duplicate_Main300_geometries':dup,'mean_difference':float(a.mean()-r.mean()),'median_difference':float(np.median(a)-np.median(r)),'exact_MW_two_sided_p':float(mannwhitneyu(a,r,method='exact',alternative='two-sided').pvalue),'exact_MW_one_sided_exploratory_p':float(mannwhitneyu(a,r,method='exact',alternative='greater').pvalue),'rank_arithmetic_trapezoidal_spearman':float(spearmanr(t.verified_J,t.trapezoidal_J).statistic),'top5_overlap':len(set(t.nsmallest(5,'arithmetic_rank').candidate_id)&set(t.nsmallest(5,'trapezoidal_rank').candidate_id)),'trapezoidal_best':t.nlargest(1,'trapezoidal_J').iloc[0].candidate_id}
(ROOT/'reviewer_analysis_summary.json').write_text(json.dumps(stats,indent=2));print(pd.DataFrame(rows).to_string(index=False));print(json.dumps(stats,indent=2))
