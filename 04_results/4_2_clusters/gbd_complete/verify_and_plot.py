"""Check fitted results / numerical DTW, and regenerate supplementary panels."""
import sys,importlib.util,itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.metrics import adjusted_rand_score
spec=importlib.util.spec_from_file_location('run',Path(__file__).with_name('run_all_models.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
out=Path(sys.argv[1]);res={}
x=np.array([0.,1.,3.]);y=np.array([0.,2.,3.]);cost=(x[:,None]-y[None,:])**2
vals=[]
def rec(i,j,v):
 if i==2 and j==2: vals.append(v);return
 for di,dj,w in [(1,1,2),(1,0,1),(0,1,1)]:
  if i+di<3 and j+dj<3:rec(i+di,j+dj,v+w*cost[i+di,j+dj])
rec(0,0,cost[0,0]);res['dtw_bruteforce_error']=abs(m.dtw(x,y)-np.sqrt(min(vals)/6))
assert res['dtw_bruteforce_error']<1e-12
f=pd.read_csv(out/'fit_diagnostics_all_models.csv');assert len(f)==126 and f.accepted.all()
res['fits_passed']=126;res['max_gradient']=float(f.max_abs_scaled_gradient.max())
res['max_multistart_eta_difference']=float(f.max_eta_difference_3_starts.max())
k=pd.read_csv(out/'kappa_all_models.csv');b=pd.read_csv(out/'parameters_all_models.csv')
mem=pd.read_csv(out/'cluster_membership.csv');assert len(mem)==126
errors=[]
for (sex,cause,model),kk in k[k.model.str.endswith('LC2')].groupby(['sex','cause','model']):
 K=kk.pivot(index='year',columns='factor',values='kappa').values
 B=b[(b.sex==sex)&(b.cause==cause)&(b.model==model)].pivot(index='age',columns='factor',values='beta').values
 th=np.pi/12;R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
 errors.append(np.max(abs(B@K.T-(B@R)@(K@R).T)))
res['rotation_surface_max_error']=float(max(errors))
(out/'validation.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
z=pd.read_csv(out/'standardized_kappa.csv');pal=plt.get_cmap('tab10')
for sex in ['Male','Female']:
 fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True,sharey=True)
 for row,model in enumerate(['CBD','Binomial LC2']):
  mm=mem[(mem.sex==sex)&(mem.model==model)]
  for col in range(3):
   ax=axes[row,col];part=mm[mm.cluster==col+1]
   if part.empty:ax.axis('off');continue
   for _,rr in part.iterrows():
    zz=z[(z.sex==sex)&(z.model==model)&(z.cause==rr.cause)]
    ax.plot(zz.year,zz.z,label=rr.code,color=pal(m.CAUSES.index(rr.cause)),lw=1.7)
   ax.set_title(f'{model} | Group {col+1}');ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.suptitle(f'{sex} | GBD Hungary | standardized period factors')
 fig.tight_layout();fig.savefig(out/f'figures7_8_{sex}.png',dpi=160);plt.close(fig)
pairs=pd.read_csv(out/'pairwise_distances.csv')
for sex in ['Male','Female']:
 fig,axes=plt.subplots(2,4,figsize=(15,8))
 for ax,model in zip(axes.flat,m.MODELS):
  D=np.zeros((9,9))
  for _,r in pairs[(pairs.sex==sex)&(pairs.model==model)].iterrows():
   i=m.CAUSES.index(r.cause1);j=m.CAUSES.index(r.cause2);D[i,j]=D[j,i]=r.dtw_p2
  ax.imshow(D,vmin=0,vmax=2,cmap='viridis')
  ax.set_xticks(range(9),m.SHORT,rotation=90);ax.set_yticks(range(9),m.SHORT);ax.set_title(model)
 axes.flat[-1].axis('off');fig.tight_layout();fig.savefig(out/f'heatmaps_{sex}.png',dpi=150);plt.close(fig)
