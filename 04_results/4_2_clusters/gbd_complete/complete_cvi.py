"""Recompute seven medoid-based cluster validity indices and compare random PAM starts."""
import sys
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.metrics import adjusted_rand_score
from run_all_models import choose,distance,CAUSES,medoids,SEED
p=Path(sys.argv[1]);z=pd.read_csv(p/'standardized_kappa.csv');mem=pd.read_csv(p/'cluster_membership.csv');old=pd.read_csv(p/'cluster_summary.csv')
rows=[];summary=[];checks=[];rng=np.random.default_rng(SEED)
for _,ss in old.iterrows():
 sex,model=ss.sex,ss.model
 Z=z[(z.sex==sex)&(z.model==model)].pivot(index='cause',columns='year',values='z').loc[CAUSES].values;d=distance(Z)
 tab,lab,_=choose(d,Z)
 before=mem[(mem.sex==sex)&(mem.model==model)].set_index('cause').loc[CAUSES].cluster.to_numpy()-1
 assert adjusted_rand_score(lab,before)==1
 best=tab.iloc[0];summary.append(dict(sex=sex,model=model,factor=ss.factor,**best.to_dict()))
 tab['sex']=sex;tab['model']=model;rows+=tab.to_dict('records')
 for k in range(2,9):
  exact=d[:,medoids(d,k)[1]].min(1).sum()
  for start in range(3):
   med=sorted(rng.choice(9,k,replace=False));obj=d[:,med].min(1).sum()
   while True:
    bestobj=obj;bestmed=med
    for a in med:
     for b in set(range(9))-set(med):
      cand=sorted((set(med)-{a})|{b});v=d[:,cand].min(1).sum()
      if v<bestobj-1e-12:bestobj=v;bestmed=cand
    if bestobj>=obj-1e-12:break
    med=bestmed;obj=bestobj
   checks.append(dict(sex=sex,model=model,K=k,start=start+1,objective=obj,global_objective=exact,gap=obj-exact))
pd.DataFrame(rows).to_csv(p/'cluster_validity_scores.csv',index=False)
pd.DataFrame(summary).to_csv(p/'cluster_summary.csv',index=False)
pd.DataFrame(checks).to_csv(p/'pam_three_start_check.csv',index=False)
print('CVI candidates:',len(rows),'PAM starts:',len(checks),'max gap:',max(x['gap'] for x in checks))
