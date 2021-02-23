import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import Counter

import allel
from ripser import ripser


class cocycleIndividualPlot:
    def __init__(self,
                 popinfo_path,
                 vcf_file=None,
                 gt_matrix_PCs=None,
                 result=None):

        self.popinfo_path = popinfo_path
        self.vcf_file = vcf_file
        self.gt_matrix_PCs = gt_matrix_PCs
        self.result = result

        if vcf_file is None and gt_matrix_PCs is None:
            raise TypeError("Either vcf_file or gt_matrix_PCs required as input")
            
    def preprocess(self):
        if self.gt_matrix_PCs is None:
            gt_matrix = self.process_vcf().get('gt_matrix')
            
            # Normalize gt_matrix by site
            gt_matrix_norm = gt_matrix - np.mean(gt_matrix,axis=1)[:, np.newaxis]
            
            # PCA
            start = time.time()
            u, s, vh = np.linalg.svd(gt_matrix_norm.T, full_matrices=False)
            print('SVD took {} secs'.format(time.time() - start))
            self.gt_matrix_PCs = np.diag(s) @ vh  
            
        # Get relevant objects from result of ripser
        if self.result is None:
            print("Getting ripser object")
            self.result = ripser(self.gt_matrix_PCs, coeff=2, maxdim=1, do_cocycles=True)
            
    def process_vcf(self):
        '''
        Reference: https://github.com/AI-sandbox/genTools/blob/943662cc9433e8943cc536435eea1ffa87127879/file_processing.py#L113
        File from Devang Agrawal to take a vcf file (plink format) and return numpy arrays
        '''
        vcf = allel.read_vcf(self.vcf_file)
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

        return {'gt_matrix': gt_matrix,
                'rs_IDs': rs_IDs,
                'ind_IDs': ind_IDs,
                'positions': positions}
                

    def display_cocycle_charts(self,
                               cocycle_number_list,
                               colordict=None,
                               cocycle_individuals_file=None,
                               birth_death_coordinates_file=None):
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
        # Obtain necessary inputs if missing
        self.preprocess()
        
        if colordict is None:
            colordict = {
                -1: 'lightgrey',
                0: 'red',
                1: 'blue',
                2: 'green',
                3: 'pink',
                4: 'black',
                5: 'yellow',
                6: 'lightblue'
            }

        diagrams = self.result['dgms']
        cocycles = self.result['cocycles']
        dgm1 = diagrams[1]
        num_cocycles = len(cocycles[1])
        num_samples = len(self.result['idx_perm'])

        # Config for plot
        xlabel, ylabel = "Birth", "Death"
        x_down, x_up = 0, np.max(dgm1)*1.2
        ax_color=np.array([0.0, 0.0, 0.0])
        point_size=0.25
        point_size_large = 1.5
        alpha = 0.8

        # Get an array of cocyles, ordered by birth death time, decreasing
        ordered_cocycle_indices = np.argsort(dgm1[:, 1] - dgm1[:, 0])[::-1]
        active_cocycle_indices = ordered_cocycle_indices[cocycle_number_list]

        # Set labels for each cocycle, -1 is for not in any active cocycles
        labels_bd = np.full(num_cocycles, -1)
        for num, idx in zip(cocycle_number_list, active_cocycle_indices):
            labels_bd[idx] = num

        # array of colors according to manually set colordict 
        colors_bd = [colordict[x] for x in labels_bd]
        sizes_bd = [point_size if x==-1 else point_size_large for x in labels_bd]

        # initialize plot
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(2, 2)
        ax2 = fig.add_subplot(gs[0, :])
        ax1 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1:, 1])

        # Get Birth Death Plot
        pd.DataFrame(dgm1).\
            plot.\
            scatter(0, 1, ax=ax1, c=colors_bd, s=sizes_bd, colorbar=False, alpha=alpha)

        ax1.set_xlim(left=x_down, right=x_up)
        ax1.set_ylim(bottom=x_down, top=x_up)
        ax1.plot([x_down, x_up], [x_down, x_up], "--", c=ax_color)
        ax1.set_aspect('equal', 'box')
        ax1.set_title('Birth Death Plot:\n location of cocycle')
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)

        # Get PC2 plot
        labels_pop = pd.read_csv(self.popinfo_path ,sep=',')
        population_count_dict = dict()
        cocycle_individuals_dict = dict()
        df_populations = pd.DataFrame(columns=['POP', 'cocycle number'])

        # Set labels for each and sample, -1 is for not in any active cocycles
        labels_pc = np.full(num_samples, -1)

        # reverse because we want label to be most persistant cocycle 
        # if individual is in multiple cocycles label is overwritten
        for i, idx in enumerate(np.flip(active_cocycle_indices)):
            cocycle = cocycles[1][idx]
            individuals = cocycle[:,:2].flatten()
            idx_flipped = cocycle_number_list[::-1][i]
            labels_pc[individuals] = idx_flipped

            # populations df for bar chart.
            # x // 2 as both haplotypes are in the popinfo but combined in the gt_matrix
            cocycle_pops = labels_pop.loc[{x // 2 for x in individuals}, 'POP']
            population_count_dict[idx_flipped] = Counter(cocycle_pops)
            cocycle_individuals_dict[idx_flipped] = cocycle_pops.index.to_list()

        # array of colors according to manually set colordict 
        colors_pc = [colordict[x] for x in labels_pc if x != -1]

        # Do background first so it appears behind
        pd.DataFrame(self.gt_matrix_PCs[labels_pc==-1]).\
            plot.\
            scatter(0, 1, ax=ax2, c=colordict[-1], s=point_size, colorbar=False, alpha=alpha)

        pd.DataFrame(self.gt_matrix_PCs[labels_pc!=-1]).\
            plot.\
            scatter(0, 1, ax=ax2, c=colors_pc, s=point_size_large, colorbar=False, alpha=alpha)
        ax2.set_title('Principal Components of cocycles')

        # get histogram of populations
        population_count_df = pd.DataFrame(population_count_dict)
        colormap = map(colordict.get, population_count_df.columns)
        population_count_df.plot.bar(stacked=True, ax=ax3, color=colormap)

        ax3.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        ax3.set_title('Counts of Populations by cocycle')
        
        # save text files
        if cocycle_individuals_file is not None:
            with open(cocycle_individuals_file, "w") as f:
                f.write(str(cocycle_individuals_dict))
                
        if birth_death_coordinates_file is not None:
            with open(birth_death_coordinates_file, "w") as f:
                np.savetxt(f, dgm1.astype(int), fmt='%i')
            

        return fig
