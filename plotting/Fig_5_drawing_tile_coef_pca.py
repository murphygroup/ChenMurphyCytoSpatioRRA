import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import seaborn as sns
import sys

if __name__ == '__main__':
    # Define parameters
    data_dir = Path(sys.argv[1])
    tissue_list = ['SPLEEN', 'THYMUS', 'LN', 'SI', 'LI']
    tile_side_list = ['2500x2500 tiles', '5000x5000 tiles', 'original images']
    tile_side_original_list = ['2500', '5000', '1']
    
    # Initialize figure and data frame
    fig, axes = plt.subplots(2, 3, figsize=(14, 10))
    hetero_df = pd.DataFrame(columns=["Heterogeneity", "Tissue", "Tile type"])
    
    # Load coefficient data for each tile type and tissue
    coef_list = []
    for tile_side, tile_side_original in zip(tile_side_list, tile_side_original_list):
        for tissue in tissue_list:
            TMC = 'Stanford' if tissue in ['LI', 'SI'] else 'Florida'
            file_path = (data_dir / TMC / tissue / 'coef_tissue.csv') if tile_side == 'original images' else \
                        (data_dir / TMC / tissue / 'random_tile' / f'coef_tile_{tile_side_original}.csv')
            coef_df = pd.read_csv(file_path)
            coef_df['tissue'] = tissue
            coef_df['side'] = tile_side
            coef_list.append(coef_df)

    coef = pd.concat(coef_list).fillna(0)

    # Standardize data and apply PCA
    scaler = StandardScaler().fit(coef.iloc[:, :-2])
    coef_pca = PCA(n_components=2).fit_transform(scaler.transform(coef.iloc[:, :-2]))
    coef[['PC1', 'PC2']] = coef_pca

    # Compute heterogeneity for each tissue and tile type
    for tile_side in tile_side_list:
        for tissue in tissue_list:
            TMC = 'Stanford' if tissue in ['LI', 'SI'] else 'Florida'
            coef_tissue = coef[(coef.tissue == tissue) & (coef.side == tile_side)].iloc[:, :-4].values
            coef_tissue_ori = pd.read_csv(data_dir / TMC / tissue / 'random_tile' / 'coef_5_100-500_1_total_across3.csv').values.reshape(1, -1)
            hetero = np.median(euclidean_distances(coef_tissue, coef_tissue_ori))
            hetero_df.loc[len(hetero_df.index)] = [hetero, tissue, tile_side]

    # Plot PCA scatter plots
    for plot_idx, tissue in enumerate(tissue_list):
        TMC = 'Stanford' if tissue in ['LI', 'SI'] else 'Florida'
        coef_tissue = coef[coef.tissue == tissue]
        ax = axes[plot_idx // 3, plot_idx % 3]
        sns.scatterplot(ax=ax, data=coef_tissue, x='PC1', y='PC2', hue='side', legend=plot_idx == 0, alpha=0.7)
        
        if plot_idx == 0:
            ax.legend(loc='upper left')
            plt.setp(ax.get_legend().get_texts(), fontsize=14)

        # Plot original image PCA point
        coef_tissue_ori = pd.read_csv(data_dir / TMC / tissue / 'random_tile' / 'coef_5_100-500_1_total_across3.csv').values.reshape(1, -1)
        coef_tissue_ori_pca = PCA(n_components=2).fit_transform(scaler.transform(coef_tissue_ori))
        ax.scatter(coef_tissue_ori_pca[:, 0], coef_tissue_ori_pca[:, 1], color='r', marker='X', alpha=0.7, label='all original images')
        
        ax.set_title(tissue, fontsize=16)
        ax.text(-31, 31, f"({chr(plot_idx + 65)})", size=16)
        ax.set_xlim(-30, 30)
        ax.set_ylim(-30, 30)
        ax.set_xlabel('PC1' if plot_idx in [2, 3, 4] else '', fontsize=16)
        ax.set_ylabel('PC2' if plot_idx in [0, 3] else '', fontsize=16)

    # Convert tile type names for heterogeneity plot
    hetero_df.replace({
        'original images': 'original\nimages',
        '5000x5000 tiles': '5000x5000\ntiles',
        '2500x2500 tiles': '2500x2500\ntiles'
    }, inplace=True)

    # Plot heterogeneity across tile types
    ax = axes[-1, -1]
    sns.lineplot(ax=ax, data=hetero_df, x='Tile type', y='Heterogeneity', hue='Tissue', palette="tab10", legend='auto')
    ax.text(-0.12, 58.5, f"({chr(plot_idx + 66)})", size=16)
    ax.set_xlabel('')
    ax.set_ylabel('Heterogeneity', fontsize=14)
    plt.setp(ax.get_legend().get_texts(), fontsize=14)

    # Save figure
    plt.savefig('./fig/tile_all_self_quad_d_no_between_dummy.png', dpi=500)
    plt.clf()
    plt.close()
