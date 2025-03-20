from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

data_dir = sys.argv[1]

def truncate_polygon(polygon, center, cutoff_distance):
    """Truncates a polygon by buffering the center point with a given cutoff distance."""
    return polygon.intersection(Point(center).buffer(cutoff_distance))

def clip_polygon_to_box(polygon, xmin, xmax, ymin, ymax):
    """Clips a polygon to the specified bounding box."""
    box = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
    return polygon.intersection(box)

# Configuration
tissue_list = ['LI', 'SI', 'THYMUS', 'LN', 'SPLEEN']
resample_percent_list = [200]
colors = ['green', 'magenta', 'cyan', 'red', 'yellow']
buffer_distance = 20
cell_num = 10000

# Define the bounding box
xmin, xmax = 800, 1800
ymin, ymax = 800, 1800

# Processing for each resample percentage and tissue type
for resample_percent in resample_percent_list:
    for tissue_idx, tissue in enumerate(tissue_list):
        TMC = 'Stanford' if tissue in ['LI', 'SI'] else 'Florida'
        file_path = f"{data_dir}/{TMC}/{tissue}/random_resample/random_simulated_pattern_5_500_1_total_across3_resample5_{resample_percent}_cell_num_{cell_num}_self_global_same_pattern_quad_d_no_between_dummy.csv"

        # Load data
        pattern_df = pd.read_csv(file_path)
        points = pattern_df.iloc[:, :2].values
        marks = pattern_df.iloc[:, 2].astype(int).values

        # Filter points within bounding box
        mask = (points[:, 0] >= xmin) & (points[:, 0] <= xmax) & (points[:, 1] >= ymin) & (points[:, 1] <= ymax)
        points, marks = points[mask], marks[mask]

        # Compute Voronoi diagram
        vor = Voronoi(points)

        # Initialize plot
        fig, ax = plt.subplots(figsize=(4, 4))
        fig.patch.set_facecolor('black')

        # Draw Voronoi cells
        for point_idx, region_idx in enumerate(vor.point_region):
            region = vor.regions[region_idx]
            if region and -1 not in region:
                polygon = [vor.vertices[i] for i in region]
                truncated = truncate_polygon(Polygon(polygon), points[point_idx], buffer_distance)
                clipped = clip_polygon_to_box(truncated, xmin, xmax, ymin, ymax)

                if not clipped.is_empty and clipped.geom_type == 'Polygon':
                    x, y = clipped.exterior.xy
                    ax.fill(x, y, color=colors[marks[point_idx]], edgecolor='white', linewidth=0.1)

        # Adjust plot limits and styling
        crop_dis = 50
        ax.set_xlim(xmin + crop_dis, xmax - crop_dis)
        ax.set_ylim(ymin + crop_dis, ymax - crop_dis)
        ax.axis('off')

        # Save figure
        save_path = f"{tissue}_global_random_resample.png"
        plt.savefig(save_path, dpi=500, bbox_inches='tight', pad_inches=0, facecolor=fig.get_facecolor(), edgecolor='none')
        plt.close()
