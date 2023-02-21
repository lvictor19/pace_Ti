import numpy as np
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler
from skcosmo.feature_selection import CUR, FPS
from skcosmo.metrics import global_reconstruction_error
from skcosmo.preprocessing import StandardFlexibleScaler
from sklearn.decomposition import TruncatedSVD

def plot_loglog_2(x,y1,y2,name,forsample=False):
    ax = plt.subplot(1,1,1)
    ax.tick_params(axis='y',which='both',direction='in')
    ax.tick_params(axis='x',which='both',direction='in')
    ax.loglog(x,y1,c='r',label='fps')
    ax.loglog(x,y2,c='b',label='cur')
    if forsample:
        ax.set_xlabel('Number of Samples Selected')
    else:
        ax.set_xlabel('Number of Features Selected')
    ax.set_ylabel(r'$R^2$')
    plt.legend()
    plt.savefig(name)
    plt.close()


def truncatedsvd_2dim(data:np.ndarray,ids:np.ndarray,plotname:str,nhighlight:int =10):
    '''
    analysis truncated svd

    Parameters
    ----------
    data : np.ndarray
    ids : np.ndarray
        sorted ids of the important features
    plotname : str
        name of the plot
    nhighlight : int  
        default : 10
        hightlight the top nhighlight points

    Returns
    -------
    None
        output plot of the data projected onto 2d space, highlighting the most important structures
    '''
    svd = TruncatedSVD(n_components=2)
    truncatedsvd=svd.fit_transform(data)
    plt.scatter(truncatedsvd[:,0],truncatedsvd[:,1])
    for _ in ids[:nhighlight]:
        plt.scatter(truncatedsvd[_,0],truncatedsvd[_,1],c='r')
    plt.savefig(plotname)
    plt.close()

def feature_selection(data:np.ndarray,plotname:str,std=False):
    '''
    input data of size (nsamples,nfeatures)
    output plot of the errors versus number of features selected

    Parameters
    ----------
    data : np.ndarray
        of size (nsamples,nfeatures)
    plotname : str
        name of the plot
    std : bool
        default : False
        using StandardScaler instead of StandardFlexibleScaler

    Returns
    -------
    None
        output plot of Loss versus number of features selected
    '''
    nfeature=len(data[0])
    if std:
        X = StandardScaler().fit_transform(data)
    else:
        X = StandardFlexibleScaler(column_wise=False).fit_transform(data)
    selector_cur = CUR(n_to_select=nfeature,progress_bar=True,score_threshold=1E-12,
                        full=False,k = 1,recompute_every = 1,tolerance=1E-12)
    selector_fps = FPS(n_to_select=nfeature,progress_bar=True,score_threshold=1E-12,
                        full=False,initialize = 0)
    id_cur=selector_cur.fit(X).selected_idx_
    id_fps=selector_fps.fit(X).selected_idx_
    errors_fps=[]
    errors_cur=[]
    k_bin=int(np.floor(nfeature/12))
    k_values=np.arange(k_bin,nfeature,k_bin)
    for k in k_values:
        errors_fps.append(global_reconstruction_error(X[:,id_fps[:k]],X))
        errors_cur.append(global_reconstruction_error(X[:,id_cur[:k]],X))
    plot_loglog_2(k_values,errors_fps,errors_cur,plotname)

def sample_selection(data:np.ndarray,plotname:str,std=False):
    '''
    input data of size (nsamples,nfeatures)
    output plot of the errors versus number of samples selected

    Parameters
    ----------
    data : np.ndarray
        of size (nsamples,nfeatures)
    plotname : str
        name of the plot
    std : bool
        default : False
        using StandardScaler instead of StandardFlexibleScaler

    Returns
    -------
    None
        output plot of Loss versus number of samples selected
        output plots of data points projected onto the plane spanned by the 
            first two important directions based on cur and fps decompositions,
    '''
    nsample=len(data)
    print(nsample)
    if std:
        X = StandardScaler().fit_transform(data.T)
    else:
        X = StandardFlexibleScaler(column_wise=False).fit_transform(data.T)
    selector_cur = CUR(n_to_select=nsample,progress_bar=True,score_threshold=1E-12,
                        full=False,k = 1,recompute_every = 1,tolerance=1E-12)
    selector_fps = FPS(n_to_select=nsample,progress_bar=True,
                        full=False,initialize = 0)
    id_cur=selector_cur.fit(X).selected_idx_
    id_fps=selector_fps.fit(X).selected_idx_
    errors_fps=[]
    errors_cur=[]
    k_bin=int(np.floor(nsample/12))
    k_values=np.arange(k_bin,nsample,k_bin)
    for k in k_values:
        errors_fps.append(global_reconstruction_error(X[:,id_fps[:k]],X))
        errors_cur.append(global_reconstruction_error(X[:,id_cur[:k]],X))
    plot_loglog_2(k_values,errors_fps,errors_cur,plotname,forsample=True)
    truncatedsvd_2dim(data,id_cur,plotname+'tsvdcur')
    truncatedsvd_2dim(data,id_fps,plotname+'tsvdfps')

