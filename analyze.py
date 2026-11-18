from pathlib import Path
import json, hashlib, platform
import numpy as np
import pandas as pd
import sklearn
from scipy.stats import mannwhitneyu
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.dummy import DummyRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.inspection import permutation_importance
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyBboxPatch
import joblib

ROOT=Path(__file__).resolve().parent
OUT=ROOT;DATA=OUT/'data';FIG=OUT/'figures'
SEED=2026
DATA.mkdir(parents=True,exist_ok=True)
FIG.mkdir(parents=True,exist_ok=True)
main=pd.read_excel(ROOT/'inputs/Paper2_ITO_Main300_ML_Phase1.xlsx',sheet_name='Merged_6000')
provided=pd.read_excel(ROOT/'inputs/Paper2_ITO_Main300_ML_Phase1.xlsx',sheet_name='Design_Splits')
cand=pd.read_excel(ROOT/'inputs/Paper2_ITO_Final_Validation_Analysis.xlsx',sheet_name='Merged_480')
original=pd.read_excel(ROOT/'inputs/Paper2_ITO_24_Candidate_COMSOL_Validation 1.xlsx',sheet_name='Candidate_Matrix_24')
reported=pd.read_excel(ROOT/'inputs/Paper2_ITO_Final_Validation_Analysis.xlsx',sheet_name='Candidate_Results')
reference=pd.read_csv(ROOT/'support/Main300_Reference_Partition_Recovery.csv')
# The user requested the original documented reference partition, not the Phase-1 reconstructed split.
mapnames={'Train':'train','Validation':'validation','Blind':'blind_test'}
main['phase1_split']=main['split'].map(mapnames)
main=main.drop(columns='split').merge(reference[['design_id','partition']],on='design_id',validate='many_to_one').rename(columns={'partition':'split'})
geometry=['ito_w_um','ito_d_um','ito_h_um','cell_w_um','cell_d_um']
features=['ito_w_um','ito_h_um','cell_w_um','cell_d_um','wavelength_um']
targets=['P_HgCdTe_W','P_ITO_W']
assert len(main)==6000 and len(cand)==480
assert main.design_id.nunique()==300 and cand.candidate_id.nunique()==24
assert not main.duplicated(['design_id','wavelength_um']).any()
assert not cand.duplicated(['candidate_id','wavelength_um']).any()
assert main.groupby('design_id').size().eq(20).all()
assert cand.groupby('candidate_id').size().eq(20).all()
for df,key in [(main,'design_id'),(cand,'candidate_id')]:
 assert np.isfinite(df[geometry+['wavelength_um','R','T','A_total']+targets]).all().all()
 assert np.allclose(df.ito_w_um,df.ito_d_um,atol=1e-12,rtol=0)
 assert df.groupby(key)[geometry].nunique().eq(1).all().all()
 assert ((df.cell_w_um-df.ito_w_um)/2>=.05-1e-12).all()
 assert ((df.cell_d_um-df.ito_d_um)/2>=.05-1e-12).all()
 grid=df.groupby(key).wavelength_um.apply(lambda s:tuple(sorted(s)))
 assert grid.nunique()==1
for g in geometry:
 left=main.groupby('design_id')[g].first().sort_index()
 right=provided.set_index('design_id')[g].sort_index()
 assert np.allclose(left,right,atol=1e-12,rtol=0)
for g in geometry:
 left=cand.groupby('candidate_id')[g].first().sort_index()
 right=original.set_index('candidate_id')[g].sort_index()
 assert np.allclose(left,right,atol=1e-12,rtol=0)
designs=main.groupby('design_id').first().reset_index()[['design_id']+geometry+['split','phase1_split']]
changed=int(designs.split.ne(designs.phase1_split).sum())
qc={'main_rows':len(main),'main_designs':300,'candidate_rows':len(cand),'candidate_designs':24,
    'split_counts':designs.split.value_counts().to_dict(),'changed_phase1_assignments':changed,
    'main_max_port':float(abs(main.R+main["T"]+main.A_total-1).max()),
    'main_max_material':float(abs(main.A_total-main.P_HgCdTe_W-main.P_ITO_W).max()),
    'candidate_max_port':float(abs(cand.R+cand["T"]+cand.A_total-1).max()),
    'candidate_max_material':float(abs(cand.A_total-cand.P_HgCdTe_W-cand.P_ITO_W).max()),
    'duplicates':0,'missing':int(main.isna().sum().sum()),'nonfinite':0,'square_violations':0,'unmatched':0}
print('QC',qc,flush=True)
main.to_csv(DATA/'Main300_Corrected_6000_Final_Split.csv',index=False)
designs.to_csv(DATA/'Main300_Final_Geometry_Partitions.csv',index=False)
cand.to_csv(DATA/'Corrected_Candidate_Spectra_480.csv',index=False)
original.to_csv(DATA/'Original_Frozen_Screening_Mean_Predictions.csv',index=False)
def metrics(y,p):
 return {'MAE':float(mean_absolute_error(y,p)),'RMSE':float(mean_squared_error(y,p)**.5),'R2':float(r2_score(y,p))}
def factory(name,seed=SEED):
 if name=='Mean':return DummyRegressor(strategy='mean')
 if name=='Polynomial Ridge':return make_pipeline(StandardScaler(),PolynomialFeatures(degree=3,include_bias=False),Ridge(alpha=1.))
 if name in ['Random Forest','Extra Trees','Wavelength-only Extra Trees']:
  cls=RandomForestRegressor if name=='Random Forest' else ExtraTreesRegressor
  return cls(n_estimators=500,max_features=1.,min_samples_leaf=1,random_state=seed,n_jobs=-1)
 if name=='Histogram Gradient Boosting':
  return HistGradientBoostingRegressor(max_iter=300,learning_rate=.08,max_leaf_nodes=31,l2_regularization=1.,random_state=seed)
 raise ValueError(name)
train=main[main.split=='train'];val=main[main.split=='validation'];blind=main[main.split=='blind_test']
assert (len(train),len(val),len(blind))==(4200,900,900)
names=['Mean','Polynomial Ridge','Random Forest','Extra Trees','Histogram Gradient Boosting','Wavelength-only Extra Trees']
results=[];selected={};models={};predictions=[];importance=[];learning=[];subsets=[]
for ti,target in enumerate(targets):
 fitted={};valid=[]
 for name in names:
  cols=['wavelength_um'] if name.startswith('Wavelength-only') else features
  model=factory(name).fit(train[cols],train[target])
  fitstats=metrics(val[target],model.predict(val[cols]))
  results.append({'target':target,'model':name,'split':'validation',**fitstats})
  fitted[name]=(model,cols)
  if not name.startswith('Wavelength-only') and name!='Mean':valid.append((fitstats['MAE'],name))
 choice=min(valid)[1];selected[target]=choice;models[target]=fitted[choice][0]
 print('SELECTED',target,choice,flush=True)
 for name,(model,cols) in fitted.items():
  bp=model.predict(blind[cols])
  results.append({'target':target,'model':name,'split':'blind',**metrics(blind[target],bp)})
  if name==choice:
   b=blind[['design_id','wavelength_um']].copy()
   b['target']=target;b['actual']=blind[target].values;b['predicted']=bp;b['signed_error']=bp-blind[target].values
   predictions.append(b)
 # Importance on validation data only, preserves blind evaluation.
 imp=permutation_importance(models[target],val[features],val[target],scoring='neg_mean_absolute_error',
                             n_repeats=30,random_state=SEED,n_jobs=1)
 for i,f in enumerate(features):
  importance.append({'target':target,'feature':f,'MAE_increase_mean':imp.importances_mean[i],'MAE_increase_sd':imp.importances_std[i]})
 ids=np.array(sorted(train.design_id.unique()))
 for n in [30,75,125,200,210]:
  for rep in range(1 if n==210 else 3):
   seed=SEED+100*n+rep
   chosen=ids if n==210 else np.sort(np.random.default_rng(seed).choice(ids,n,replace=False))
   part=train[train.design_id.isin(chosen)]
   m=factory(choice,SEED).fit(part[features],part[target])
   learning.append({'target':target,'model':choice,'training_designs':n,'replicate':rep,'subset_seed':seed,
                    'model_seed':SEED,**metrics(val[target],m.predict(val[features]))})
   for design_id in chosen:subsets.append({'target':target,'training_designs':n,'replicate':rep,'subset_seed':seed,'design_id':design_id})
 joblib.dump(models[target],OUT/(target+'_audit_model.joblib'))
 print('FINISHED',target,flush=True)
bench=pd.DataFrame(results);bp=pd.concat(predictions,ignore_index=True)
imp=pd.DataFrame(importance);lc=pd.DataFrame(learning)
bench.to_csv(DATA/'Final_Model_Benchmarks.csv',index=False);bp.to_csv(DATA/'Final_Blind_Predictions.csv',index=False)
imp.to_csv(DATA/'Validation_Permutation_Importance.csv',index=False);lc.to_csv(DATA/'Learning_Curve_Runs.csv',index=False)
pd.DataFrame(subsets).to_csv(DATA/'Learning_Curve_Subset_IDs.csv',index=False)
# Geometry-level blind objectives, using the final audit refits.
wide=bp.pivot(index=['design_id','wavelength_um'],columns='target',values=['actual','predicted'])
per=wide.groupby('design_id').mean()
per.columns=['_'.join(x) for x in per.columns];per=per.reset_index()
per['actual_J']=per['actual_P_HgCdTe_W']-per['actual_P_ITO_W'];per['predicted_J']=per['predicted_P_HgCdTe_W']-per['predicted_P_ITO_W']
per['J_error']=per.predicted_J-per.actual_J
per.to_csv(DATA/'Blind_Geometry_Objectives.csv',index=False)
verified=cand.groupby('candidate_id').agg(mean_HgCdTe=('P_HgCdTe_W','mean'),mean_ITO=('P_ITO_W','mean'),
             mean_A_total=('A_total','mean'),max_port=('energy_closure',lambda s:abs(s).max()),
             max_material=('material_balance',lambda s:abs(s).max())).reset_index()
verified['verified_J']=verified.mean_HgCdTe-verified.mean_ITO
verified=verified.merge(original,on='candidate_id',validate='one_to_one')
for name,target in [('HgCdTe','P_HgCdTe_W'),('ITO','P_ITO_W')]:
 pred=models[target].predict(cand[features])
 cand['audit_pred_'+name]=pred
 z=cand.groupby('candidate_id')['audit_pred_'+name].mean()
 verified['audit_mean_'+name]=verified.candidate_id.map(z)
 verified['original_'+name+'_signed_error']=verified['pred_mean_'+name+'_W']-verified['mean_'+name]
 verified['original_'+name+'_abs_error']=abs(verified['original_'+name+'_signed_error'])
verified['original_J_signed_error']=verified.pred_objective_J-verified.verified_J
verified['original_J_abs_error']=abs(verified.original_J_signed_error)
verified['audit_J']=verified.audit_mean_HgCdTe-verified.audit_mean_ITO
for key in ['verified_J']:
 chk=reported.set_index('candidate_id')[key]
 assert np.allclose(verified.set_index('candidate_id')[key].sort_index(),chk.sort_index(),atol=1e-12,rtol=0)
peaks=[]
for key,g in cand.groupby('candidate_id'):
 actual=g.loc[g.P_HgCdTe_W.idxmax()]
 audit=g.loc[g.audit_pred_HgCdTe.idxmax()]
 peaks.append({'candidate_id':key,'verified_peak_W':actual.P_HgCdTe_W,'verified_peak_um':actual.wavelength_um,
               'audit_pred_peak_W':audit.audit_pred_HgCdTe,'audit_pred_peak_um':audit.wavelength_um})
verified=verified.merge(pd.DataFrame(peaks),on='candidate_id',validate='one_to_one')
verified.to_csv(DATA/'Corrected_All_24_Candidate_Results.csv',index=False)
cand.to_csv(DATA/'Candidate_Spectra_With_Audit_Predictions.csv',index=False)
groups=verified.groupby(verified.candidate_id.str[:2])
guided=verified[verified.candidate_id.str.startswith('BO')].sort_values('candidate_id')
random=verified[verified.candidate_id.str.startswith('RS')].sort_values('candidate_id')
a=guided.verified_J.to_numpy();b=random.verified_J.to_numpy()
rng=np.random.default_rng(SEED)
boot=a[rng.integers(0,12,(100000,12))].mean(1)-b[rng.integers(0,12,(100000,12))].mean(1)
test=mannwhitneyu(a,b,alternative='greater',method='asymptotic',use_continuity=True)
pooled=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
best=verified.sort_values('verified_J',ascending=False).iloc[0]
best_random=random.sort_values('verified_J',ascending=False).iloc[0].candidate_id
mainavg=main.groupby('design_id')[targets].mean();mainavg['J']=mainavg.P_HgCdTe_W-mainavg.P_ITO_W
stats={'guided_mean':a.mean(),'random_mean':b.mean(),'difference':a.mean()-b.mean(),
       'relative_percent':(a.mean()-b.mean())/b.mean()*100,'CI':np.quantile(boot,[.025,.975]).tolist(),
       'archived_CI':[0.01299867151779658,0.01896927587049323],
       'U':test.statistic,'p':test.pvalue,'cohen_d':(a.mean()-b.mean())/pooled,
       'best':best.candidate_id,'best_J':best.verified_J,'best_random':best_random,
       'best_main':mainavg.J.idxmax(),'best_main_J':mainavg.J.max(),
       'original_candidate_MAE':{g:{q:float(d['original_'+q+'_abs_error'].mean()) for q in ['HgCdTe','ITO','J']} for g,d in [('guided',guided),('random',random)]},
       'blind_J':metrics(per.actual_J,per.predicted_J)}
mainavg.to_csv(DATA/'Main300_Design_Mean_Powers.csv')
summary={'qc':qc,'selected_models':selected,'statistics':stats,'seed':SEED,
 'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__},
 'features':features,'models':{t:models[t].get_params() for t in targets},
 'sources':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (ROOT/'inputs').glob('Paper2_ITO*')}}
(OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2,default=lambda v: float(v) if isinstance(v,np.generic) else repr(v)))
print(json.dumps(summary,indent=2,default=lambda v: float(v) if isinstance(v,np.generic) else repr(v)),flush=True)
