import numpy as np
import allel

import plotly
import plotly.graph_objects as go



def process_vcf(vcf_file):
    '''
    Reference: https://github.com/AI-sandbox/genTools/blob/943662cc9433e8943cc536435eea1ffa87127879/file_processing.py#L113
    File from Devang Agrawal to take a vcf file (plink format) and return numpy arrays
    '''
    vcf = allel.read_vcf(vcf_file)
    # Genotype array
    gt = vcf['calldata/GT']
    n_variants, n_samples, ploidy = gt.shape
    gt_matrix = gt.reshape(n_variants, n_samples * ploidy).astype(np.float32)
    np.place(gt_matrix, gt_matrix < 0, np.nan)
    
    # ID
    IDs = vcf['variants/ID']
    rs_IDs = [int(x.split(':')[-1]) for x in IDs]
    
    # Samples
    samples = vcf['samples']
    ind_IDs = []
    for sample in samples:
        ind_IDs.append(sample + '_A')
        ind_IDs.append(sample + '_B')
    ind_IDs = np.array(ind_IDs)
    
    # Positions
    positions = vcf['variants/POS'].tolist()
    
    return gt_matrix, rs_IDs, ind_IDs, positions


def plot_cocycle(X, cocycle, **kwargs):
    
    """
    Function from Brad to plot cocyles from Ripser
    """
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=X[:,0], y=X[:,1],
        mode='markers',
    ))
    for k in range(cocycle.shape[0]):
        [i, j, val] = cocycle[k, :]
        Xi = X[i]
        Xj = X[j]
    
        fig.add_shape(
            # Line Diagonal
                type="line",
                x0=Xi[0],
                y0=Xi[1],
                x1=Xj[0],
                y1=Xj[1],
                line=dict(
                    width=4,
                )
        )
    fig.update_layout(**kwargs)
    return fig  



def plot_cocycle_2D(X, cocycle, D, thresh, **kwargs):
    """
    plot edges in cocycle below given threshold
    X: N x 2 numpy array of point locations
    cocycle: ripser cocycle
    D: N x N distance matrix
    thresh: threshold parameter
    kwargs: passed onto figure layout
    """
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=X[:,0], y=X[:,1],
        mode='markers',
    ))
    edge_x = []
    edge_y = []
    N = X.shape[0]
    for i in range(N):
        for j in range(N):
            if D[i, j] <= thresh:
                edge_x.extend([X[i,0], X[j,0], None])
                edge_y.extend([X[i,1], X[j,1], None])
                
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')
     )

    edge_x = []
    edge_y = []
    for k in range(cocycle.shape[0]):
        [i, j, val] = cocycle[k, :]
        if D[i, j] <= thresh:
            edge_x.extend([X[i,0], X[j,0], None])
            edge_y.extend([X[i,1], X[j,1], None])
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=2, color='red'),
        hoverinfo='none',
        mode='lines')
     )
    fig.update_layout(**kwargs)
    return fig


def plot_representative_2D(X, F, R, pair, D, thresh, **kwargs):
    """
    plot representative
    X: 2-dimensional locations of points
    F: bats FilteredSimplicialComplex 
    R: bats ReducedFilteredChainComplex
    pair: bats PersistencePari
    D: N x N distance matrix
    thresh: threshold parameter
    kwargs: passed onto figure layout
    """
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=X[:,0], y=X[:,1],
        mode='markers',
    ))
    edge_x = []
    edge_y = []
    N = X.shape[0]
    for i in range(N):
        for j in range(N):
            if D[i, j] <= thresh:
                edge_x.extend([X[i,0], X[j,0], None])
                edge_y.extend([X[i,1], X[j,1], None])
                
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')
     )

    edge_x = []
    edge_y = []
    r = R.representative(pair)
    nzind = r.nzinds()
    cpx = F.complex()
    for k in nzind:
        [i, j] = cpx.get_simplex(1, k)
        if D[i, j] <= thresh:
            edge_x.extend([X[i,0], X[j,0], None])
            edge_y.extend([X[i,1], X[j,1], None])
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=2, color='red'),
        hoverinfo='none',
        mode='lines')
     )
    fig.update_layout(**kwargs)
    return fig


def sample_keeping_pop_breadth(labels, n):
    '''
    Take in the popinfo.csv file and return a list of indices representing a sample of points
    which maintains the breadth in population. It does this by taking a single row, without replacement
    from each population until there are n indices.
    '''
    result = []

    df = labels.copy()

    while True:
        if df['POP'].value_counts().count() < n:
            n -= df['POP'].value_counts().count()
            # Add one of each
            for pop in df['POP'].unique():
                idx = df[df['POP']==pop].sample(1).index
                result.extend(idx)
                df = df.drop(idx)
        else:
            # sample from remainig ones
            result.extend(df.sample(n).index)
            break
            
    return result