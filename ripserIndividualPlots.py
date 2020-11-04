import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import Counter

def display_cocycle_charts(cocycle_number_list,
                           gt_matrix_PCs, 
                           colordict,
                           result=None,
                           popinfo_path='~/../projects/HimalGenAsia/HimalGen.popinfo.csv'):
    """
    Plots population information by cocycle in the Rips complex.
    
    For a list of cocycles, plots the individuals which belong to them on PC2 and 
    Birth-Death scatterplots, and as a histogram by population frequency.
    
    Parameters
    ------------
    cocycle_number_list (list): Odered list of cocycles, zero indexed with 0 being the most 
                                persistant, ie longest birth death time.
    gt_matrix_PCs (numpy array):
    colordict(dict):            A dictionary mapping the members of cocycle_number_list
                                to colors which should contrast each other as much as possible
                                as these are the colors the cocycles will show up with in the 
                                plots. -1 is the key of the base color.
    result: (dict):             Output of ripser.ripser()
    popinfo_path (str):         The location of the popinfo csv file.
                                
  
    Returns: 
    A matplotlib fig object containing 3 subplots
    """
    # Get relevant objects from result of ripser
    if result is None:
        print("Getting ripser object")
        result = ripser(gt_matrix_PCs, coeff=2, maxdim=1, do_cocycles=True)
        
    diagrams = result['dgms']
    cocycles = result['cocycles']
    dgm1 = diagrams[1]
    num_cocycles = len(cocycles[1])
    num_samples = len(result['idx_perm'])

    # Config for plot
    xlabel, ylabel = "Birth", "Death"
    x_down, x_up = 0, np.max(dgm1)*1.2
    ax_color=np.array([0.0, 0.0, 0.0])
    point_size=0.5
    point_size_large = 5

    # Get an array of cocyles, ordered by birth death time, decreasing
    ordered_cocycle_indices = np.argsort(dgm1[:, 1] - dgm1[:, 0])[::-1]
    active_cocycle_indices = ordered_cocycle_indices[cocycle_number_list]

    # Set labels for each cocycle, -1 is for not in any active cocycles
    labels_bd = np.full(num_cocycles, -1)
    for num, idx in zip(cocycle_number_list, active_cocycle_indices):
        labels_bd[idx] = num

    # array of colors according to manually set colordict 
    colors_bd = [colordict[x] for x in labels_bd]
    sizes_bd = [point_size if x==-1 else point_size_large for x in labels_bd ]

    # initialize plot
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax2 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1:, 1])

    # Get Birth Death Plot
    pd.DataFrame(dgm1).\
        plot.\
        scatter(0,1,ax=ax1, c=colors_bd, s=sizes_bd, colorbar=False)

    ax1.set_xlim(left=x_down, right=x_up)
    ax1.set_ylim(bottom=x_down, top=x_up)
    ax1.plot([x_down, x_up], [x_down, x_up], "--", c=ax_color)
    ax1.set_aspect('equal', 'box')
    ax1.set_title('BD Plot showing location of cocycle')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    # Get PC2 plot
    labels_pop = pd.read_csv(popinfo_path ,sep=',')
    population_count_dict = dict()
    df_populations = pd.DataFrame(columns=['POP', 'cocycle number'])
    
    # Set labels for each and sample, -1 is for not in any active cocycles
    labels_pc = np.full(num_samples, -1)

    # reverse because we want label to be most persistant cocycle 
    # if individual is in multiple cocycles label is overwritten
    for i, idx in enumerate(np.flip(active_cocycle_indices)):
        cocycle = cocycles[1][idx]
        individuals = cocycle[:,:2].flatten()
        idx_flipped = len(cocycle_number_list) - i - 1
        labels_pc[individuals] = idx_flipped

        # populations df for bar chart. 
        # x // 2 as both haplotypes are in the popinfo but combined in the gt_matrix
        cocycle_pops = labels_pop.loc[{x // 2 for x in individuals}, 'POP']
        population_count_dict[idx_flipped] = Counter(cocycle_pops)

    # array of colors according to manually set colordict 
    colors_pc = [colordict[x] for x in labels_pc]
    sizes_pc = [point_size if x==-1 else point_size_large for x in labels_pc]

    pd.DataFrame(gt_matrix_PCs).\
        plot.\
        scatter(0,1, ax=ax2, c=colors_pc, s=sizes_pc, colorbar=False)
    ax2.set_title('Principal Components of cocycle')

    # get histogram of populations
    population_count_df = pd.DataFrame(population_count_dict)
    colormap = map(colordict.get, population_count_df.columns)
    population_count_df.plot.bar(stacked=True, ax=ax3, color=colormap)
    ax3.set_title('Counts of Populations in this cocycle')
    
    return fig
