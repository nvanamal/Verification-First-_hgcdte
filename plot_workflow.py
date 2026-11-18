from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
R=Path(__file__).resolve().parent
fig,ax=plt.subplots(figsize=(10,3.8));ax.set_xlim(0,10);ax.set_ylim(0,3.7);ax.axis('off')
boxes=[(.15,2.2,'Documented ITO table\n300 geometries × 20 wavelengths'),(3.55,2.2,'Geometry-grouped evaluation\nRetrospective held-out audit'),(6.95,2.2,'Square-footprint screening\nFrozen candidate sets'),(6.95,.2,'Independent COMSOL verification\n12 guided + 12 random'),(3.55,.2,'Targeted numerical refinement\nM201 settings recorded separately'),(.15,.2,'Performance interpretation\nObjective-dependent group separation')]
for x,y,t in boxes:
 ax.add_patch(FancyBboxPatch((x,y),2.85,1.15,boxstyle='round,pad=.03',ec='#267078',fc='#edf3f3',lw=1.2));ax.text(x+1.425,y+.575,t,ha='center',va='center',fontsize=10)
for a,b in [((3.07,2.775),(3.47,2.775)),((6.47,2.775),(6.87,2.775)),((8.375,2.12),(8.375,1.43)),((6.87,.775),(6.47,.775)),((3.47,.775),(3.07,.775))]:ax.add_patch(FancyArrowPatch(a,b,arrowstyle='->',mutation_scale=14,color='#267078',lw=1.5))
fig.tight_layout(pad=.3)
for ext in ['png','pdf']:fig.savefig(R/'figures'/('workflow.'+ext),dpi=240,bbox_inches='tight')
