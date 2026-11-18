from pathlib import Path
import io,json,hashlib,shutil
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
RAW=ROOT/'refinement_exports';RAW.mkdir(exist_ok=True)
DATA=ROOT/'data';DATA.mkdir(exist_ok=True)
# Only archived refinement_exports are used in the public package.
# Explicit mapping: the original bo06_*.csv exports contain 480 production rows, despite their model filenames.
MAPPING={'BO06_finemesh_port.csv':'20 wave port.csv','BO06_finemesh_hgcdte.csv':'20 wave hgcdte.csv','BO06_finemesh_ito.csv':'20 wave ito(1).csv','BO06_40wl_port.csv':'40 wave port.csv','BO06_40wl_hgcdte.csv':'40 wave hgcdte.csv','BO06_40wl_ito.csv':'40 wave ito.csv','M201_finemesh_port.csv':'m201 port.csv','M201_finemesh_hgcdte.csv':'m201hgcdte.csv','M201_finemesh_ito.csv':'m201 ito.csv'}
def read(path):
 return pd.read_csv(io.StringIO('\n'.join(l for l in path.read_text().splitlines() if not l.startswith('%'))),header=None)
manifest=[]
for dst,src in MAPPING.items():
 path=RAW/dst
 manifest.append({'archive_name':dst,'original_filename':src,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(DATA/'Refinement_Source_Manifest.csv',index=False)
rows=[]
for name in ['BO06_finemesh','BO06_40wl','M201_finemesh']:
 h,i,p=[read(RAW/(name+'_'+kind+'.csv')) for kind in ['hgcdte','ito','port']]
 assert np.allclose(h[0],i[0]) and np.allclose(h[0],p[0]);assert not h[0].duplicated().any()
 df=pd.DataFrame({'wavelength_um':h[0],'P_HgCdTe_W':h[1],'P_ITO_W':i[1],'R':p[1],'T':p[2],'A_total':p[4]});df['J_W']=df.P_HgCdTe_W-df.P_ITO_W
 df.to_csv(DATA/(name+'_merged.csv'),index=False)
 j=int(df.P_HgCdTe_W.idxmax());rows.append({'record':name,'wavelengths':len(df),'mean_HgCdTe_W':df.P_HgCdTe_W.mean(),'mean_ITO_W':df.P_ITO_W.mean(),'arithmetic_J_W':df.J_W.mean(),'trapezoidal_J_W':np.trapezoid(df.J_W,df.wavelength_um)/(df.wavelength_um.max()-df.wavelength_um.min()),'peak_HgCdTe_W':df.loc[j,'P_HgCdTe_W'],'peak_wavelength_um':df.loc[j,'wavelength_um'],'max_material_residual':abs(df.A_total-df.P_HgCdTe_W-df.P_ITO_W).max(),'max_port_residual':abs(df.R+df['T']+df.A_total-1).max()})
s=pd.DataFrame(rows);s.to_csv(DATA/'Refinement_Summary.csv',index=False)
f=s.iloc[0];q=s.iloc[1]
summary={'same_fine_mesh_spectral_change_percent':float(100*(q.trapezoidal_J_W-f.trapezoidal_J_W)/f.trapezoidal_J_W),'record_notes':'The 20-point BO06 files are fine-mesh baseline exports; bo06_*.csv contain the 480-row production candidate set. The prior 0.6931% mixed production-mesh and fine-mesh records.','ito_filename':'ito.txt','ito_sha256':json.loads((ROOT/'Material_Provenance.json').read_text())['sha256'],'ito_wavelength_unit':'micrometer numerical values, supplied through author-confirmed 1[1/um]*c_const/freq'}
(ROOT/'refinement_summary.json').write_text(json.dumps(summary,indent=2));print(s.to_string(index=False));print(summary)
