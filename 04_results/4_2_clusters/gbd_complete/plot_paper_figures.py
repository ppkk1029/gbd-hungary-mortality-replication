"""Regenerate individual GBD analogues of the paper's Figures 5/6/7 from results/ CSV."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
P=Path(__file__).resolve().parent/'results';P.mkdir(exist_ok=True)
z=pd.read_csv(P/'standardized_kappa.csv');mem=pd.read_csv(P/'cluster_membership.csv')
dist=pd.read_csv(P/'pairwise_distances.csv')
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.15})
def series(sex,model,cause):
 d=z[(z.sex==sex)&(z.model==model)&(z.cause==cause)].sort_values('year')
 assert len(d)==36
 return d.year.to_numpy(),d.z.to_numpy()
for sex,n,other,code in [('Male',5,'Digestive diseases','DIG'),('Female',6,'Diabetes and kidney diseases','DKD')]:
 base='Cardiovascular diseases';x,a=series(sex,'Poisson LC2',base);_,b=series(sex,'Poisson LC2',other)
 path=pd.read_csv(P/f'alignment_path_{sex}.csv')
 i=(path.year_1-1980).to_numpy();j=(path.year_2-1980).to_numpy()
 d=dist[(dist.sex==sex)&(dist.model=='Poisson LC2')&(((dist.cause1==base)&(dist.cause2==other))|((dist.cause2==base)&(dist.cause1==other)))]
 assert len(d)==1
 dtw=float(d.dtw_p2.iloc[0]);local=np.abs(a[:,None]-b[None,:])
 fig=plt.figure(figsize=(14,7));gs=GridSpec(2,3,figure=fig,width_ratios=[1.4,1,1],height_ratios=[1.1,.8])
 ax=fig.add_subplot(gs[:,0]);ax.plot(x,a,label='CVD / kappa(2)');ax.plot(x,b,label=code+' / kappa(2)')
 ax.legend();ax.set(title='Two standardized period factors',xlabel='Year',ylabel='z-score')
 ax=fig.add_subplot(gs[0,1:]);im=ax.imshow(local,origin='lower',aspect='auto',extent=[1979.5,2015.5,1979.5,2015.5])
 ax.plot(x[j],x[i],color='white',lw=1.6)
 ax.set(title=f'Local absolute difference / path; DTW={dtw:.3f}',xlabel=code+' year',ylabel='CVD year')
 fig.colorbar(im,ax=ax,fraction=.04)
 ax=fig.add_subplot(gs[1,1:]);ax.plot(x,a,label='CVD');ax.plot(x,b+4,label=code+' + 4')
 for u,v in zip(i[::2],j[::2]):
  ax.plot([x[u],x[v]],[a[u],b[v]+4],color='gray',alpha=.3,lw=.6)
 ax.legend();ax.set(title='DTW correspondences',xlabel='Year')
 fig.suptitle(f'Our Figure {n}: GBD Hungary / {sex} / Poisson LC2')
 fig.tight_layout();fig.savefig(P/f'Figure_{n}_GBD_{sex}_Poisson_LC2.png',dpi=160);plt.close(fig)
m=mem[(mem.sex=='Male')&(mem.model=='CBD')]
fig,axes=plt.subplots(1,2,figsize=(14,5.2),sharey=True)
for k in [1,2]:
 ax=axes[k-1]
 for _,r in m[m.cluster==k].iterrows():
  x,y=series('Male','CBD',r.cause);ax.plot(x,y,label=r.code,lw=1.8)
 ax.set(title=f'Cluster {k}',xlabel='Year',ylabel='Standardized CBD first period factor')
 ax.legend(ncol=2)
fig.suptitle('Our Figure 7: GBD Hungary / Male / CBD kappa(1)')
fig.tight_layout();fig.savefig(P/'Figure_7_GBD_Male_CBD_Clusters.png',dpi=160);plt.close(fig)

# Original Figure 8: male two-factor Binomial LC second period factor.
m=mem[(mem.sex=='Male')&(mem.model=='Binomial LC2')]
fig,axes=plt.subplots(1,2,figsize=(14,5.2),sharey=True)
for k in [1,2]:
 ax=axes[k-1]
 for _,rr in m[m.cluster==k].iterrows():
  x,y=series('Male','Binomial LC2',rr.cause);ax.plot(x,y,label=rr.code,lw=1.8)
 ax.set(title=f'Cluster {k}',xlabel='Year',ylabel='Standardized Binomial LC2 second period factor')
 ax.legend(ncol=2);ax.grid(alpha=.15)
fig.suptitle('Our Figure 8: GBD Hungary / Male / Binomial Lee–Carter kappa(2)')
fig.tight_layout();fig.savefig(P/'Figure_8_GBD_Male_Binomial_LC2.png',dpi=160);plt.close(fig)
