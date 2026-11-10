"""GBD Hungary | seven single-population models + Section 4.2 DTW clustering.
Reproducibility source recovered from actual 2026-11-10 full-results archive.
Usage: OPENBLAS_NUM_THREADS=1 python run_all_models.py fitted_rates.csv results
The supplied CSV must contain ORIGINAL deaths_gbd and exposure_proxy; fitted_count is not input.
Dependencies: numpy pandas scipy matplotlib scikit-learn. Python 3.10+.
"""
import sys, json, hashlib, itertools
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, xlogy
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score, adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
SEED=20261110
MODELS=['Poisson LC1','Binomial LC1','Poisson LC2','Binomial LC2','CBD','Plat2','Plat3']
CAUSES=['Cardiovascular diseases','Neoplasms','Digestive diseases','Chronic respiratory diseases','Neurological disorders','Unintentional injuries','Diabetes and kidney diseases','Self-harm and interpersonal violence','Respiratory infections and tuberculosis']
SHORT=['CVD','NEO','DIG','CRD','NEU','INJ','DKD','SHV','RIT']
AGES=np.arange(62,88,5); YEARS=np.arange(1980,2016)

def canonical(a,B,K):
    center=K.mean(0); a=a+B@center; K=K-center
    u,s,vt=np.linalg.svd(B@K.T,full_matrices=False);r=B.shape[1]
    B=u[:,:r];K=vt[:r].T*s[:r]
    for j in range(r):
        sign=1 if B[np.argmax(abs(B[:,j])),j]>=0 else -1
        B[:,j]*=sign;K[:,j]*=sign
    return a,B,K,s

def fit(D,E,model):
    A,T=D.shape;binom=model.startswith('Binomial') or model=='CBD'
    N=E+D/2 if binom else E
    z=np.log(D/(N-D)) if binom else np.log(D/E)
    total=D.sum();rng=np.random.default_rng(SEED)
    lc='LC' in model
    if lc:
        r=int(model[-1]);a0=z.mean(1);u,s,vt=np.linalg.svd(z-a0[:,None],full_matrices=False)
        combos=list(itertools.combinations(range(A),r));anchor=max(combos,key=lambda ix:abs(np.linalg.det(u[list(ix),:r])))
        M=u[list(anchor),:r];B0=u[:,:r]@np.linalg.inv(M);K0=(vt[:r].T*s[:r])@M.T
        free=[i for i in range(A) if i not in anchor];nb=len(free)*r
        def pack(a,B,K):return np.r_[a,B[free].ravel(),K[:-1].ravel()]
        def unpack(v):
            a=v[:A];B=np.zeros((A,r));B[list(anchor)]=np.eye(r);B[free]=v[A:A+nb].reshape(-1,r)
            K=np.zeros((T,r));K[:-1]=v[A+nb:].reshape(-1,r);K[-1]=-K[:-1].sum(0)
            return a,B,K
        start=pack(a0,B0,K0)
    else:
        r=3 if model=='Plat3' else 2
        B0=np.c_[np.ones(A),((AGES.mean()-AGES)/10)]
        if model=='CBD':B0[:,1]*=-1
        if r==3:B0=np.c_[B0,np.maximum((AGES.mean()-AGES)/10,0)]
        a0=np.zeros(A) if model=='CBD' else z.mean(1)
        K0=np.linalg.lstsq(B0,z-a0[:,None],rcond=None)[0].T
        def pack(a,B,K):return K.ravel() if model=='CBD' else np.r_[a,K[:-1].ravel()]
        def unpack(v):
            if model=='CBD':return np.zeros(A),B0,v.reshape(T,r)
            K=np.zeros((T,r));K[:-1]=v[A:].reshape(-1,r);K[-1]=-K[:-1].sum(0)
            return v[:A],B0,K
        start=pack(a0,B0,K0)
    def calc(v):
        a,B,K=unpack(v);eta=a[:,None]+B@K.T
        if binom:mu=N*expit(eta);val=np.sum(N*np.logaddexp(0,eta)-D*eta)/total
        else:mu=E*np.exp(np.clip(eta,-700,50));val=np.sum(mu-D*eta)/total
        res=(mu-D)/total;ga=res.sum(1);gk=res.T@B
        if lc:g=np.r_[ga,(res@K)[free].ravel(),(gk[:-1]-gk[-1]).ravel()]
        elif model=='CBD':g=gk.ravel()
        else:g=np.r_[ga,(gk[:-1]-gk[-1]).ravel()]
        return val,g
    starts=[start]+[start+rng.normal(0,.015,len(start)) for _ in range(2)]
    fits=[minimize(calc,v,jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':4000}) for v in starts]
    best=min(fits,key=lambda f:f.fun);a,B,K=unpack(best.x);eta=a[:,None]+B@K.T
    grad=float(max(abs(calc(best.x)[1])))
    if lc:a,B,K,s=canonical(a,B,K)
    else:s=np.linalg.svd(B@K.T,compute_uv=False)
    mu=N*expit(eta) if binom else E*np.exp(eta)
    dev=2*np.sum(xlogy(D,D/mu)-D+mu) if not binom else 2*np.sum(xlogy(D,D/mu)+xlogy(N-D,(N-D)/(N-mu)))
    rate=expit(eta)/(1-expit(eta)/2) if binom else np.exp(eta)
    diag={'deviance':dev,'max_abs_scaled_gradient':grad,'accepted':grad<1e-6,'objective_range_3_starts':max(f.fun for f in fits)-min(f.fun for f in fits),'max_eta_difference_3_starts':max(np.max(abs((unpack(f.x)[0][:,None]+unpack(f.x)[1]@unpack(f.x)[2].T)-eta)) for f in fits),'mape_percent':100*np.mean(abs(rate-D/E)/(D/E)),'factor2_factor1_ratio':s[1]/s[0] if len(s)>1 else np.nan,'reconstruction_error':np.max(abs(a[:,None]+B@K.T-eta))}
    assert diag['accepted'],(model,diag)
    return a,B,K,diag

def dtw(x,y,p=2,window=None,path=False):
    n=len(x);m=len(y);cost=np.abs(x[:,None]-y[None,:])**p
    acc=np.full((n,m),np.inf);prev=np.full((n,m),-1,int);acc[0,0]=cost[0,0]
    for i in range(n):
        for j in range(m):
            if i==0 and j==0:continue
            if window is not None and abs(i-j)>window:continue
            opts=[acc[i-1,j-1]+2*cost[i,j] if i and j else np.inf,acc[i-1,j]+cost[i,j] if i else np.inf,acc[i,j-1]+cost[i,j] if j else np.inf]
            k=np.argmin(opts);acc[i,j]=opts[k];prev[i,j]=k
    d=(acc[-1,-1]/(n+m))**(1/p)
    if not path:return d
    pts=[];i=n-1;j=m-1
    while True:
        pts.append((i,j))
        if i==0 and j==0:break
        k=prev[i,j];i-=k in (0,1);j-=k in (0,2)
    return d,cost,acc,np.array(pts[::-1])

def distance(Z,p=2,window=None):
    n=len(Z);d=np.zeros((n,n))
    for i in range(n):
        for j in range(i):d[i,j]=d[j,i]=dtw(Z[i],Z[j],p,window)
    return d

def medoids(d,k):
    # n=9 makes exact enumeration feasible, avoiding local PAM optima.
    med=min(itertools.combinations(range(len(d)),k),key=lambda ix:d[:,ix].min(1).sum())
    lab=d[:,med].argmin(1)
    for j,i in enumerate(med):lab[i]=j
    return lab,np.array(med)

def metrics(d,lab,Z):
    groups=[np.where(lab==j)[0] for j in np.unique(lab)];n=len(d);k=len(groups)
    meds=[g[np.argmin(d[np.ix_(g,g)].sum(1))] for g in groups]
    diam=max(d[np.ix_(g,g)].max() for g in groups)
    sep=min(d[np.ix_(a,b)].min() for i,a in enumerate(groups) for b in groups[:i])
    scatter=np.array([d[g,m].mean() for g,m in zip(groups,meds)])
    mm=d[np.ix_(meds,meds)];rat=(scatter[:,None]+scatter[None,:])/np.maximum(mm,1e-15);np.fill_diagonal(rat,-np.inf)
    globalmed=np.argmin(d.sum(1));within=sum(np.sum(d[g,m]**2) for g,m in zip(groups,meds))
    between=sum(len(g)*d[m,globalmed]**2 for g,m in zip(groups,meds))
    centers=np.array([Z[g].mean(0) for g in groups]);grand=Z.mean(0)
    we=sum(np.sum((Z[g]-c)**2) for g,c in zip(groups,centers));be=sum(len(g)*np.sum((c-grand)**2) for g,c in zip(groups,centers))
    cop=sum(d[g,m].sum()/max(np.min(np.max(d[np.ix_(np.setdiff1d(np.arange(n),g),g)],axis=0)),1e-15) for g,m in zip(groups,meds))/n
    dbstar=np.mean([max(scatter[i]+np.delete(scatter,i))/max(np.min(np.delete(mm[i],i)),1e-15) for i in range(k)])
    linear_between=sum(len(g)*d[m,globalmed] for g,m in zip(groups,meds));linear_within=sum(d[g,m].sum() for g,m in zip(groups,meds))
    sf=-np.expm1(-np.exp(np.clip(linear_between/(n*k)-scatter.sum(),-700,6)))
    return {'silhouette':silhouette_score(d,lab,metric='precomputed'),'dunn':sep/max(diam,1e-15),'DB_medoid':np.max(rat,axis=1).mean(),'DBstar_medoid':dbstar,'CH_dtwclust_medoid':(n-k)/(k-1)*linear_between/linear_within,'SF_medoid':sf,'CH_squared_medoid':(between/(k-1))/(within/(n-k)),'CH_euclidean':(be/(k-1))/(we/(n-k)),'COP_medoid':cop,'objective':linear_within}

def choose(d,Z):
    rows=[];labs={};L=linkage(squareform(d),method='average')
    for k in range(2,len(d)):
        for method in ['average','kmedoids_exact']:
            lab=cut_tree(L,n_clusters=k).ravel() if method=='average' else medoids(d,k)[0]
            rows.append({'method':method,'K':k,**metrics(d,lab,Z)});labs[(method,k)]=lab
    tab=pd.DataFrame(rows).sort_values(['silhouette','K','method'],ascending=[False,True,True])
    best=tab.iloc[0];return tab,labs[(best['method'],int(best.K))],L

def main():
    if len(sys.argv)!=3:raise SystemExit('Usage: python run_all_models.py fitted_rates.csv output_dir')
    inp=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
    raw=pd.read_csv(inp);assert len(raw)==3888 and not raw.duplicated(['sex_name','cause_name','age_name','year']).any()
    assert set(raw.year)==set(YEARS) and set(raw.cause_name)==set(CAUSES)
    assert (raw.deaths_gbd>0).all() and (raw.exposure_proxy>0).all()
    kappas=[];params=[];diagnostics=[];allfits={}
    for sex in ['Male','Female']:
      for cause in CAUSES:
        f=raw[(raw.sex_name==sex)&(raw.cause_name==cause)].sort_values(['age_name','year'])
        D=f.deaths_gbd.to_numpy().reshape(6,36);E=f.exposure_proxy.to_numpy().reshape(6,36)
        for model in MODELS:
            a,B,K,d=fit(D,E,model);allfits[(sex,cause,model)]=(a,B,K)
            diagnostics.append(dict(sex=sex,cause=cause,model=model,**d))
            for j in range(K.shape[1]):
              for t,y in enumerate(YEARS):kappas.append(dict(sex=sex,cause=cause,model=model,factor=j+1,year=y,kappa=K[t,j]))
              for i,x in enumerate(AGES):params.append(dict(sex=sex,cause=cause,model=model,factor=j+1,age=x,alpha=a[i],beta=B[i,j]))
        print('fit',sex,cause,flush=True)
    pd.DataFrame(diagnostics).to_csv(out/'fit_diagnostics_all_models.csv',index=False)
    pd.DataFrame(kappas).to_csv(out/'kappa_all_models.csv',index=False)
    pd.DataFrame(params).to_csv(out/'parameters_all_models.csv',index=False)
    summary=[];membership=[];scores=[];sens=[];pairs=[];zrows=[]
    for sex in ['Male','Female']:
      fig,axes=plt.subplots(4,2,figsize=(14,16));fig2,axes2=plt.subplots(4,2,figsize=(14,15))
      for mi,model in enumerate(MODELS):
        factor=1 if model.endswith('LC2') else 0
        X=np.array([allfits[(sex,c,model)][2][:,factor] for c in CAUSES])
        Z=(X-X.mean(1,keepdims=True))/X.std(1,ddof=1,keepdims=True)
        d=distance(Z);tab,lab,L=choose(d,Z);best=tab.iloc[0]
        summary.append(dict(sex=sex,model=model,factor=factor+1,**best.to_dict()))
        for i,c in enumerate(CAUSES):
            membership.append(dict(sex=sex,model=model,cause=c,code=SHORT[i],cluster=int(lab[i])+1))
            for t,y in enumerate(YEARS):zrows.append(dict(sex=sex,model=model,cause=c,year=y,z=Z[i,t]))
            for j in range(i):pairs.append(dict(sex=sex,model=model,cause1=c,cause2=CAUSES[j],dtw_p2=d[i,j]))
        tab['sex']=sex;tab['model']=model;scores+=tab.to_dict('records')
        ax=axes.flat[mi];dendrogram(L,labels=SHORT,ax=ax,color_threshold=0)
        ax.set_title(f'{model}: K={int(best.K)}, silhouette={best.silhouette:.3f}')
        ax2=axes2.flat[mi]
        for i in range(9):ax2.plot(YEARS,Z[i],label=f'{SHORT[i]} (G{lab[i]+1})',lw=1.5,alpha=.85)
        ax2.set_title(model);ax2.legend(ncol=3,fontsize=6);ax2.set_ylabel('standardized kappa')
        for tag,p,w in [('window5',2,5),('L1_legacy',1,None)]:
            ds=distance(Z,p,w);ts,ls,_=choose(ds,Z)
            sens.append(dict(sex=sex,model=model,variant=tag,K=int(ts.iloc[0].K),silhouette=ts.iloc[0].silhouette,ARI=adjusted_rand_score(lab,ls)))
        if sex=='Female':
            keep=[i for i in range(9) if i!=6]
            ts,ls,_=choose(d[np.ix_(keep,keep)],Z[keep])
            sens.append(dict(sex=sex,model=model,variant='without_DKD',K=int(ts.iloc[0].K),silhouette=ts.iloc[0].silhouette,ARI=adjusted_rand_score(lab[keep],ls)))
        if model.endswith('LC2'):
            for angle in [15,30,45]:
                th=np.deg2rad(angle);R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
                XR=np.array([(allfits[(sex,c,model)][2]@R)[:,1] for c in CAUSES])
                ZR=(XR-XR.mean(1,keepdims=True))/XR.std(1,ddof=1,keepdims=True)
                ts,ls,_=choose(distance(ZR),ZR)
                sens.append(dict(sex=sex,model=model,variant=f'rotation_{angle}',K=int(ts.iloc[0].K),silhouette=ts.iloc[0].silhouette,ARI=adjusted_rand_score(lab,ls)))
        if model=='Poisson LC2':
            i,j=(0,2) if sex=='Male' else (0,6)
            dist,cost,acc,pts=dtw(Z[i],Z[j],path=True)
            fp,ap=plt.subplots(1,3,figsize=(15,4))
            ap[0].plot(YEARS,Z[i],label=SHORT[i]);ap[0].plot(YEARS,Z[j],label=SHORT[j]);ap[0].legend();ap[0].set_title(f'{sex}: Poisson LC2 / kappa2')
            ap[1].imshow(np.sqrt(cost),origin='lower',extent=[1980,2015,1980,2015],aspect='auto')
            ap[1].plot(YEARS[pts[:,1]],YEARS[pts[:,0]],color='white')
            ap[1].set_title(f'Local absolute difference / path; DTW={dist:.3f}')
            ap[2].plot(YEARS,Z[i],label=SHORT[i]);ap[2].plot(YEARS,Z[j]+4,label=SHORT[j]+' + 4')
            for u,v in pts[::2]:ap[2].plot([YEARS[u],YEARS[v]],[Z[i,u],Z[j,v]+4],color='gray',alpha=.3,lw=.6)
            ap[2].set_title('Optimal alignment')
            fp.tight_layout();fp.savefig(out/f'alignment_{sex}.png',dpi=160);plt.close(fp)
            pd.DataFrame({'year_1':YEARS[pts[:,0]],'year_2':YEARS[pts[:,1]]}).to_csv(out/f'alignment_path_{sex}.csv',index=False)
      for ax in [axes.flat[-1],axes2.flat[-1]]:ax.axis('off')
      fig.suptitle(f'GBD Hungary | {sex} | average-linkage dendrograms')
      fig.tight_layout();fig.savefig(out/f'dendrograms_{sex}.png',dpi=150);plt.close(fig)
      fig2.suptitle(f'GBD Hungary | {sex} | selected clusters')
      fig2.tight_layout();fig2.savefig(out/f'clusters_{sex}.png',dpi=150);plt.close(fig2)
    for name,rows in [('cluster_summary',summary),('cluster_membership',membership),('cluster_validity_scores',scores),('sensitivity',sens),('pairwise_distances',pairs),('standardized_kappa',zrows)]:
        pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False)
    (out/'manifest.json').write_text(json.dumps({'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'seed':SEED,'years':[1980,2015],'ages':AGES.tolist(),'models':MODELS,'n_fits':126,'dtw':'symmetric2 weighted squared local differences; sqrt(cost/(n+m)); z-score ddof=1','selection':'max silhouette across K=2..8 and average/exhaustive k-medoids; tie lower K','binomial_exposure':'E0=Ec+D_c/2 approximation','lc_identification':'SVD orthogonal loadings, descending singular values, largest absolute loading positive'},indent=2))
    print(pd.DataFrame(summary)[['sex','model','K','silhouette']].to_string(index=False))
if __name__=='__main__':main()
