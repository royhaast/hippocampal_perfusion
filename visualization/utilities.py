def downsample_data(data, sampling_factor=0.5, anti_aliasing=True, round=True):
    import numpy as np
    from skimage.transform import rescale

    flatmap_dims = [254,126]

    data_unfolded = data.reshape(
        flatmap_dims[1],flatmap_dims[0]
    )

    data_unfolded_rescaled = rescale(
        data_unfolded, sampling_factor, anti_aliasing=anti_aliasing, multichannel=False, mode='reflect'
    )

    flatmap_dims_new = [data_unfolded_rescaled.shape[0], data_unfolded_rescaled.shape[1]]
    data_rescaled = data_unfolded_rescaled.reshape(
        int(flatmap_dims_new[1])*int(flatmap_dims_new[0]))

    if round:
        data_rescaled = np.round(data_rescaled)

    return data_rescaled

def subfield_stats(all_data, i, idx, subjects, subfields, within='Subfield'):
    import numpy as np
    import pandas as pd
    import pingouin as pg
    from brainspace.utils.parcellation import reduce_by_labels    
    
    # Prepare data for stats
    data = np.hstack((
        np.vstack((
            np.array([ [s]*5 for s in subjects]).flatten(),
            ['Left']*5*len(subjects),
            ['Sub','CA1','CA2','CA3','CA4/DG']*len(subjects),
            reduce_by_labels(all_data[:,:,0,i], subfields, axis=1).flatten(order='F')
        )),
        np.vstack((
            np.array([ [s]*5 for s in subjects]).flatten(),
            ['Right']*5*len(subjects),
            ['Sub','CA1','CA2','CA3','CA4/DG']*len(subjects),
            reduce_by_labels(all_data[:,:,1,i], subfields, axis=1).flatten(order='F')
        ))
    ))

    df_stats = pd.DataFrame(
        data=data.T,
        columns=['Subject','Hemisphere','Subfield','Data']
    )
    df_stats = df_stats.astype({'Data': float}) 

    # Non-parametric repeated measures ANOVA
    rm = pg.friedman(
        dv='Data', within=within, subject='Subject', data=df_stats)
    
    # FDR-corrected pairwise comparisons
    pw = pg.pairwise_tests(
        data=df_stats, dv='Data', within=within, parametric=False,
        padjust='fdr_bh', subject='Subject'
    )

    # # Collect significant pairs
    # pairs   = []
    # pvalues = []
    # for j in range(0,len(pw)):
    #     if (pw.iloc[j]['p-corr'] <= .05):
    #         pairs.append((pw.iloc[j]['A'],pw.iloc[j]['B']))
    #         pvalues.append(pw.iloc[j]['p-corr'])
            
    return df_stats, rm, pw

def plot_subfield_data(metric_data, subfields, metric_info, fname,
                       scale=1, stats=False, stats_data=[],
                       stats_results=[], save_fig=False):
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from statannotations.Annotator import Annotator
    from brainspace.utils.parcellation import reduce_by_labels
    
    fwidth  = 5.*scale # total width of the figure in inches
    fheight = 3. # total height of the figure in inches

    fig = plt.figure(figsize=(fwidth, fheight))

    left_margin   = 0.95 / fwidth
    right_margin  = 0.2 / fwidth
    bottom_margin = 0.5 / fheight
    top_margin    = 0.25 / fheight

    x = left_margin    # horiz. position of bottom-left corner
    y = bottom_margin  # vert. position of bottom-left corner
    w = 1 - (left_margin + right_margin) # width of axes
    h = 1 - (bottom_margin + top_margin) # height of axes

    ax = fig.add_axes([x, y, w, h])

    xloc =  0.25 / fwidth 
    yloc =  y + h / 2.  
     
    ax.set_ylabel(f'{metric_info[3]} {metric_info[4]}', fontsize='12', fontweight='bold',
                 verticalalignment='top', horizontalalignment='center')    
    ax.yaxis.set_label_coords(xloc, yloc, transform = fig.transFigure)

    for h, hemi in enumerate(['Left','Right']):
        data = np.mean(metric_data[:,:,h], axis=1)

        # Subfield average per vertex across subjects and hemispheres
        x_data = downsample_data(subfields, 0.5, False, True)
        y_data = downsample_data(data, 0.5, True, False)

        sns.stripplot(x=x_data, y=y_data, alpha=.05, jitter=0.25, ax=ax)

        # Subfield average per subject and hemisphere
        x_data = ['Sub','CA1','CA2','CA3','CA4/DG']*metric_data.shape[1]
        y_data = reduce_by_labels(metric_data[:,:,h], subfields, axis=1).flatten(order='F')
        sns.stripplot(
            x=x_data,y=y_data, marker='o' if h == 0 else 'D', 
            size=7, linewidth=1, edgecolor='black', zorder=10, ax=ax
        ) 
    
    # Display average lines
    ax.set_ylim(metric_info[5])
    bottom, top = ax.get_ylim()
    avg_padding = (top-bottom)*.05        
    for s, ss in enumerate(np.unique(subfields)):
        subfield_mean = np.nanmean(metric_data[subfields==ss,:,:])
        subfield_std = np.nanstd(metric_data[subfields==ss,:,:])
        xmin  = 0.05+(0.20*s)
        xmax  = 0.20+(0.20*s)
        xtext = 0.3+s

        # # Create and add S.D. box
        # data_to_figure = ax.transData.transform
        # data_to_axes = lambda x: ax.transAxes.transform(data_to_figure(x))
        # rect = patches.Rectangle((50, 100), 40, 30, linewidth=0, facecolor='black', alpha=.2)
        # ax.add_patch(rect)

        ax.axhline(y=subfield_mean, xmin=xmin, xmax=xmax, color='black', clip_on=False, zorder=5)
        ax.annotate('%.2f' % subfield_mean, (xtext, subfield_mean + avg_padding), xycoords='data',
                   fontsize=12, fontstyle='italic', rotation=90, zorder=10)

    # Add stats
    if stats:
        annot = Annotator(ax=ax, pairs=stats_results[0], data=stats_data, x='Subfield', y='Data')
        annot.configure(
            test=None, test_short_name='pw', loc='inside', line_height=0,
            text_offset=-2, line_offset_to_group=0)
        annot.set_pvalues(pvalues=stats_results[1])
        annot.annotate()
        
    ax.set_xticklabels(['Sub','CA1','CA2','CA3','CA4/\nDG'], fontsize='12')
    ax.tick_params(axis='y', labelsize=12)
    sns.despine()

    if save_fig:
        fig.savefig(fname, dpi=600, bbox_inches='tight', transparent=True)
    plt.show(block=False)

def lineplot(data, subfields, colors, ylim, ylabel, fname):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from brainspace.utils.parcellation import reduce_by_labels
    
    fwidth  = 3.5 #2.75 # total width of the figure in inches
    fheight = 3  # total height of the figure in inches

    fig = plt.figure(figsize=(fwidth, fheight))

    # Margins
    left_margin   = 0.95 / fwidth
    right_margin  = 0.2 / fwidth
    bottom_margin = 0.5 / fheight
    top_margin    = 0.25 / fheight

    x = left_margin    # horiz. position of bottom-left corner
    y = bottom_margin  # vert. position of bottom-left corner
    w = 1 - (left_margin + right_margin) # width of axes
    h = 1 - (bottom_margin + top_margin) # height of axes

    ax = fig.add_axes([x, y, w, h])

    xloc =  0.25 / fwidth 
    yloc =  y + h / 2.  
    ax.set_ylabel(ylabel, fontsize='12', fontweight='bold',
                 verticalalignment='top', horizontalalignment='center')    
    ax.yaxis.set_label_coords(xloc, yloc, transform = fig.transFigure)

    # Subfield average per subject and hemisphere
    x_data = ['Sub','CA1','CA2','CA3','CA4/DG']*data.shape[1]
    for h in [0,1]:
        if len(colors) > 1:
            for c, cc in zip(range(0,data.shape[1]),colors):
                y_data = reduce_by_labels(data[:,:,h,c], subfields, axis=1).flatten()
                sns.lineplot(x=x_data, y=y_data, color=cc, linestyle='-' if h == 0 else '--', ax=ax)
        else:
            y_data = reduce_by_labels(data[:,:,h], subfields, axis=1).flatten()
            sns.lineplot(x=x_data, y=y_data, color=colors[0], linestyle='-' if h == 0 else '--', ax=ax)        

    ax.set_xticklabels(['Sub','CA1','CA2','CA3','CA4/\nDG'], fontsize='12')
    ax.set_ylim(ylim)
    ax.tick_params(axis='y', labelsize=12)
    sns.despine()
    
    plt.show(block=False)
    fig.savefig(fname, dpi=600, bbox_inches='tight', transparent=True)
    
def scatterplot(xdata, ydata, xlabel, ylabel, xlim, ylim, fname):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats import linregress, spearmanr
    from sklearn.preprocessing import MinMaxScaler
    
    fwidth  = 2.75 
    fheight = 3

    fig = plt.figure(figsize=(fwidth, fheight))

    # Margins
    left_margin   = 0.8 / fwidth
    right_margin  = 0.2 / fwidth
    bottom_margin = 0.5 / fheight
    top_margin    = 0.25 / fheight

    # Create axes
    x = left_margin   
    y = bottom_margin
    w = 1 - (left_margin + right_margin)
    h = 1 - (bottom_margin + top_margin)

    ax = fig.add_axes([x, y, w, h])

    # Define positions
    xloc =  0.25 / fwidth 
    yloc =  y + h / 2.  

    ydata = ydata.T
    xlabel = xlabel
    ylabel = ylabel

    min_max_scaler = MinMaxScaler()
    annotation = ''
    if ylabel == 'PVE':
        for i, (tissue, color) in enumerate(zip(['GM','WM','CSF'],['red','green','blue'])):
            
            y_scaled = min_max_scaler.fit_transform(ydata[i,:][:,np.newaxis])

            sns.regplot(
                y=y_scaled, x=xdata,
                scatter_kws={'ec': 'None', 's': 1, 'color': color, 'alpha': 0.1},
                line_kws={'color': color},
                ax=ax
            )

            r, p = spearmanr(xdata, ydata[i,:])
            pperm = correlation_between_maps(xdata, ydata[i,:], 10000)
            annotation = annotation + (
                "{0}: R$^2$ = {1} ({2})\n".format(tissue, "%.2f" % r**2, pval_denoter(pperm, prefix=True))
            )
            
        # Decorate Y axis
        ax.set_yticks([0,1])
        ax.set_yticklabels(['Min','Max'])            
    else:
        sns.regplot(
            y=ydata, x=xdata,
            scatter_kws={'ec': 'None', 's': 1, 'color': 'black', 'alpha': 0.1},
            line_kws={'color': 'black'},
            ax=ax
        )

        r, p = spearmanr(xdata, ydata)
        pperm = correlation_between_maps(xdata, ydata, 10000)
        annotation = annotation + (
            "R$^2$ = {0} ({1})\n".format("%.2f" % r**2, pval_denoter(pperm, prefix=True))
        )    

    t = ax.annotate(annotation[:-1], xy=(1,1), xycoords='axes fraction', ha='right', va='top',
                    color='black', bbox=dict(facecolor='white', alpha=0.5, linewidth=0))

    # Decorate plot
    ax.set_xlabel(xlabel, fontsize='12', fontweight='bold')
    ax.set_ylabel(ylabel, fontsize='12', fontweight='bold',
                 verticalalignment='top', horizontalalignment='center')    
    ax.yaxis.set_label_coords(xloc, yloc, transform = fig.transFigure)
    ax.tick_params(axis='both', labelsize=12)
    sns.despine()

    plt.show(block=False)
    fig.savefig(fname, dpi=600, bbox_inches='tight', transparent=True)
    
def plot_subfield_pairs(df, rm, pw, i, metric_info, fname, scale=1, save_fig=False):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # Relabel
    pval_bin = pd.cut(
        pw['p-corr'],
        [0, .001, .005, .01, .05, 1],
        labels=['∗∗∗∗', '∗∗∗', '∗∗', '∗', '']
    )

    # Calculate averages per subfield (across subjects and hemispheres)
    subfield_averages = df.groupby(['Subfield'])['Data'].mean()

    # Prepare data for plotting as heatmap
    annot_matrix = np.repeat('', 5*5).reshape(5,5).astype(object)
    stats_matrix = np.zeros((5,5))

    ind_dict = {'Sub': 0, 'CA1': 1, 'CA2': 2, 'CA3': 3, 'CA4/DG': 4}

    for j in range(0,len(pw)):
        idx  = ind_dict[pw['A'][j]]
        idy  = ind_dict[pw['B'][j]]
        annot_matrix[idx,idy] = pval_bin.values[j]
        stats_matrix[idx,idy] = subfield_averages[pw['A'][j]]-subfield_averages[pw['B'][j]]

    stats_matrix[1:,1:] = np.rot90(np.fliplr(stats_matrix[1:,1:]))*-1
    annot_matrix[1:,1:] = np.rot90(np.fliplr(annot_matrix[1:,1:]))
    stats_matrix = stats_matrix[1:,:]
    annot_matrix = annot_matrix[1:,:]
    idx = np.nonzero(stats_matrix[:,:])[0]
    idy = np.nonzero(stats_matrix[:,:])[1]    

    # Plot
    # fig, ax = plt.subplots(1,1, figsize=(5,3))

    fwidth  = 5.*scale # total width of the figure in inches
    fheight = 3. # total height of the figure in inches

    fig = plt.figure(figsize=(fwidth, fheight))

    left_margin   = 0.95 / fwidth
    right_margin  = 0.2 / fwidth
    bottom_margin = 0.5 / fheight
    top_margin    = 0.25 / fheight

    x = left_margin    # horiz. position of bottom-left corner
    y = bottom_margin  # vert. position of bottom-left corner
    w = 1 - (left_margin + right_margin) # width of axes
    h = 1 - (bottom_margin + top_margin) # height of axes

    ax = fig.add_axes([x, y, w, h])

    xloc =  0.25 / fwidth 
    yloc =  y + h / 2.  
    
    vmin = metric_info[6][0]
    vmax = metric_info[6][1]

    hm = sns.heatmap(
        stats_matrix[:,:],
        annot=annot_matrix[:,:], fmt='',
        annot_kws={'fontsize': 10},
        center=0, cmap='bwr',
        vmin=vmin, vmax=vmax,
        linewidth=2,
        cbar_kws={
            "orientation": "horizontal",
            "shrink": .84, "fraction": .1, "pad": .3
        },
        ax=ax
    )

    hm.set_yticklabels(['CA1','CA2','CA3','CA4/\nDG'], fontsize='12', rotation=0) 
    hm.set_xticklabels(['Sub','CA1','CA2','CA3','CA4/\nDG'], fontsize='12', rotation=0) 
    # hm.get_xaxis().set_visible(False)

    cbar = ax.collections[0].colorbar
    cbar.set_ticks([vmin,vmax])
    cbar.ax.tick_params(labelsize=10, length=0)
    cbar.set_label(
        f"Difference {metric_info[4]}", weight='bold', fontsize=10, labelpad=-10)

    if save_fig:
        fig.savefig(fname, dpi=600, bbox_inches='tight')
    fig.show()

def generate_permsamples(nperm):
    import numpy as np
    
    np.random.seed(1234)
    
    # Generate distribution of shifts in both directions
    rand_x = np.random.randint(0, high=126, size=nperm)
    rand_y = np.random.randint(0, high=254, size=nperm)
    
    return rand_x, rand_y

def correlation_between_maps_permutation(out, x, y, rot, shift_pd, shift_ap, p):
    from scipy.ndimage import rotate, shift
    from netneurotools.stats import efficient_pearsonr
    
    # Rotate map
    x_rot = rotate(
        x, rot, axes=(1, 0), reshape=False, output=None,
        order=3, mode='wrap', cval=0.0, prefilter=True
    )

    # Translate map
    x_rot_shift = shift(
        x_rot, [shift_pd, shift_ap], output=None, order=3, 
        mode='wrap', cval=0.0, prefilter=True
    )

    out[p] = efficient_pearsonr(x_rot_shift.flatten(), y)[0]    

def correlation_between_maps(x, y, nperm, return_perms=False):
    import os
    import numpy as np
    from scipy.ndimage import rotate, shift
    from joblib import Parallel, delayed, dump, load
    from netneurotools.stats import efficient_pearsonr

    # To 2D plane
    x_2d = x.reshape((126,254), order='C')

    # Generate distribution of shifts in both directions    
    rotation = np.random.randint(1,360,nperm)
    shift_pd = np.random.randint(-63,64,nperm)
    shift_ap = np.random.randint(-127,128,nperm)
    
    # Empirical Pearson's correlation
    emp_corr = efficient_pearsonr(x, y)[0]
    
    # Obtain permutated Pearson's correlations
    # Temporary location for dumping data
    tmp_dir = os.environ['SLURM_TMPDIR']

    # Dump output to file for memory-efficiency
    output_filename_memmap = os.path.join(tmp_dir, f'output_memmap')
    output_memmap = np.memmap(
        output_filename_memmap, dtype=float, mode='w+',
        shape=(nperm)
    )

    # Fill data array
    Parallel(n_jobs=-1)(
        delayed(correlation_between_maps_permutation)(
            output_memmap, x_2d, y, rotation[p], shift_pd[p], shift_ap[p], p
        ) for p in range(nperm)
    )   

    # Save to Numpy file
    perm_corr = np.array(output_memmap)    

    # Calculate permuted p-value
    pval = np.mean(np.abs(perm_corr) >= np.abs(emp_corr))          
    
    if return_perms:
        return emp_corr, pval, perm_corr
    else:
        return emp_corr, pval

def pval_denoter(pval, prefix=False):
    if .01 < pval <= .05:
        return '{} < .05'.format('p' if prefix else '')
    elif .005 < pval <= .01:
        return '{} < .01'.format('p' if prefix else '')
    elif .001 < pval <= .005:
        return '{} < .005'.format('p' if prefix else '')
    elif pval <= .001:
        return '{} < .001'.format('p' if prefix else '')
    else:
        return 'n.s.'    
    
def plot_data_test(xdata, ydata, subfields, xdata_info, ydata_info, fname):
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from statannotations.Annotator import Annotator
    from brainspace.utils.parcellation import reduce_by_labels
    
    fwidth  = 2.75 # total width of the figure in inches
    fheight = 3. # total height of the figure in inches

    fig = plt.figure(figsize=(fwidth, fheight))

    # Margins
    left_margin   = 0.8 / fwidth
    right_margin  = 0.2 / fwidth
    bottom_margin = 0.5 / fheight
    top_margin    = 0.25 / fheight

    x = left_margin    # horiz. position of bottom-left corner
    y = bottom_margin  # vert. position of bottom-left corner
    w = 1 - (left_margin + right_margin) # width of axes
    h = 1 - (bottom_margin + top_margin) # height of axes

    ax = fig.add_axes([x, y, w, h])

    xloc =  0.25 / fwidth 
    yloc =  y + h / 2.  

    ax.set_xlabel(f'{xdata_info[3]} {xdata_info[4]}', fontsize='12', fontweight='bold',
                 verticalalignment='top', horizontalalignment='center')   
    ax.set_ylabel(f'{ydata_info[3]} {ydata_info[4]}', fontsize='12', fontweight='bold',
                 verticalalignment='top', horizontalalignment='center')    
    ax.yaxis.set_label_coords(xloc, yloc, transform = fig.transFigure)

    x_concat = []
    y_concat = []
    for h, hemi in enumerate(['Left','Right']):
        # Subfield average per subject and hemisphere
        x_data = reduce_by_labels(xdata[:,:,h], subfields, axis=1).flatten(order='F')
        y_data = reduce_by_labels(ydata[:,:,h], subfields, axis=1).flatten(order='F')
        hue_data = ['Sub','CA1','CA2','CA3','CA4/DG']*xdata.shape[1]
        sns.scatterplot(
            x=x_data,y=y_data, hue=hue_data, legend=False, ax=ax
        ) 
        
        x_concat = np.concatenate((x_concat, x_data))
        y_concat = np.concatenate((y_concat, y_data))
    
    sns.regplot(
        x=x_concat,y=y_concat, scatter=False, line_kws={'color': 'black'}, ax=ax
    ) 
        
    ax.tick_params(axis='y', labelsize=12)
    sns.despine()

    fig.savefig(fname, dpi=600, bbox_inches='tight', transparent=True)
    plt.show(block=False)    
    
def bootstrap_iteration_full(data_in, data_out, i, idx, idy, subfields, run_index):
    import numpy as np
    from utilities import subfield_stats
    
    for r in range(data_out.shape[1]):
        run_mask = np.argwhere(np.isin(run_index, idx[:r+1])).reshape(-1)

        for s in range(data_out.shape[0]):
            vertex_data = data_in[
                np.ix_(
                    np.arange(0,32004), idy[:s+1], np.array([0,1]), run_mask
                )
            ]
            
            # Average
            avg = np.nanmean(vertex_data, axis=(3,2,1))
            data_out[s,r,:,0,i] = avg
            
            # Standard variation
            std = np.nanstd(np.nanmean(vertex_data, axis=(3,2)), axis=1)
            data_out[s,r,:,1,i] = std

            # tSNR
            tsnr = np.nanmean(np.nanmean(vertex_data, axis=3)/np.nanstd(vertex_data, axis=3), axis=(2,1))
            data_out[s,r,:,2,i] = tsnr

            # Calculate between subfield statistic
            vertex_data = np.nanmean(vertex_data, axis=-1)[:,:,:,np.newaxis]

            _, rm, _ = subfield_stats(vertex_data, 0, 1, idy[:s+1], subfields)
            data_out[s,r,0,3,i] = rm['Q'].values[0]
            data_out[s,r,1,3,i] = rm['p-unc'].values[0]
            
    return data_out