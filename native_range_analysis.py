"""Native-ITO-range sensitivity from archived production-grid power spectra.
Linearly interpolate power at native table bounds and use normalized trapezoids.
This does not rerun the electromagnetic model or remove extrapolation influence.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent
lo,hi=2.623526293,4.712033503
s=pd.read_csv(R/'data/Corrected_Candidate_Spectra_480.csv');rows=[]
for cid,g in s.groupby('candidate_id',sort=True):
 g=g.sort_values('wavelength_um');x=g.wavelength_um.to_numpy();q=np.r_[lo,x[(x>lo)&(x<hi)],hi]
 vals={}
 for col in ['P_HgCdTe_W','P_ITO_W']:
  y=g[col].to_numpy();vals[col]=float(np.trapezoid(np.interp(q,x,y),q)/(hi-lo))
 j=float((g.P_HgCdTe_W-g.P_ITO_W).mean());full=float(np.trapezoid(g.P_HgCdTe_W-g.P_ITO_W,x)/(x[-1]-x[0]))
 rows.append(dict(candidate_id=cid,candidate_type=g.candidate_type.iloc[0],native_mean_HgCdTe_W=vals['P_HgCdTe_W'],native_mean_ITO_W=vals['P_ITO_W'],native_J_W=vals['P_HgCdTe_W']-vals['P_ITO_W'],full_arithmetic_J_W=j,full_trapezoidal_J_W=full))
d=pd.DataFrame(rows)
for name in ['native_J_W','full_arithmetic_J_W','full_trapezoidal_J_W']:d[name+'_rank']=d[name].rank(ascending=False,method='min').astype(int)
d.to_csv(R/'data/Native_ITO_Range_Sensitivity.csv',index=False)
a=d[d.candidate_id.str.startswith('BO')];b=d[d.candidate_id.str.startswith('RS')]
r=dict(bounds_um=[lo,hi],method='Piecewise-linear spectral interpolation at bounds, normalized trapezoidal integration; existing production-grid responses only',best_candidate=d.loc[d.native_J_W.idxmax(),'candidate_id'],guided_mean=float(a.native_J_W.mean()),random_mean=float(b.native_J_W.mean()),difference=float(a.native_J_W.mean()-b.native_J_W.mean()),spearman_vs_full_arithmetic=float(spearmanr(d.native_J_W,d.full_arithmetic_J_W).statistic),spearman_vs_full_trapezoidal=float(spearmanr(d.native_J_W,d.full_trapezoidal_J_W).statistic),all_144_pairs_favor_guided=bool(a.native_J_W.min()>b.native_J_W.max()),top_five_overlap=len(set(d.nlargest(5,'native_J_W').candidate_id)&set(d.nlargest(5,'full_arithmetic_J_W').candidate_id)))
(R/'native_range_summary.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
