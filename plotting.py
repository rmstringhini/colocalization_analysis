import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from skimage import io, filters, morphology
from skimage.util import img_as_float
import seaborn as sns
from pathlib import Path
from matplotlib.colors import LogNorm
from matplotlib.colors import LogNorm, PowerNorm
from matplotlib.colors import LinearSegmentedColormap



def plot_cytofluorogram(self, a, b, ch1_thresh=None, ch2_thresh=None,
                     save_path=None, bins=256):
    plt.style.use('dark_background')
    ch1 = self.ch1_raw.flatten()
    ch2 = self.ch2_raw.flatten()

    # Remove zero-zero pixels 
    valid = (ch1 + ch2) > 0
    ch1 = ch1[valid]
    ch2 = ch2[valid]

    fig, ax = plt.subplots(figsize=(7, 7))

    h = ax.hist2d(
        ch1,
        ch2,
        bins=bins,
        norm=LogNorm(),
        cmap="hot"
    )

  #  plt.colorbar(h[3], ax=ax, label="Pixel count (log)")

    # Regression line
    x_vals = np.array([0, ch1.max()])
    y_vals = a * x_vals + b
    ax.plot(x_vals, y_vals, color='cyan', lw=2)


    ax.set_xlabel("Channel 1 Intensity")
    ax.set_ylabel("Channel 2 Intensity")
    ax.set_title("Cytofluorogram (Coloc2 style)")
    #ax.legend()
    ax.grid(False)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()
    
def plot_analysis(self, save_path=None, save_individual=True):
        fig = plt.figure(figsize=(16, 12))

        if save_individual and save_path:
            output_dir = Path(save_path).parent / (Path(save_path).stem + "_individual")
            output_dir.mkdir(exist_ok=True)

        ax1 = plt.subplot(3, 4, 1)
        ax1.imshow(self.ch1_raw, cmap='Reds')
        ax1.set_title('Channel 1', fontsize=12, fontweight='bold')
        ax1.axis('off')

        ax2 = plt.subplot(3, 4, 2)
        ax2.imshow(self.ch2_raw, cmap='Greens')
        ax2.set_title('Channel 2', fontsize=12, fontweight='bold')
        ax2.axis('off')

        ax3 = plt.subplot(3, 4, 3)
        merged = np.zeros((*self.ch1_raw.shape, 3))
        merged[..., 0] = self.ch1_raw  # Red
        merged[..., 1] = self.ch2_raw  # Green
        ax3.imshow(merged)
        ax3.set_title('Merged', fontsize=12, fontweight='bold')
        ax3.axis('off')

        ax4 = plt.subplot(3, 4, 4)
        coloc_mask = self.ch1_thresh & self.ch2_thresh
        ax4.imshow(coloc_mask, cmap='hot')
        ax4.set_title('Colocalized Pixels', fontsize=12, fontweight='bold')
        ax4.axis('off')

        ax5 = plt.subplot(3, 4, 5)
        ch1_flat = self.ch1_raw[self.mask]
        ch2_flat = self.ch2_raw[self.mask]
        # ch1_flat = self.ch1_raw.flatten()
        # ch2_flat = self.ch2_raw.flatten()

        if len(ch1_flat) > 10000:
            idx = np.random.choice(len(ch1_flat), 10000, replace=False)
            ch1_plot = ch1_flat[idx]
            ch2_plot = ch2_flat[idx]
        else:
            ch1_plot = ch1_flat
            ch2_plot = ch2_flat

        ax5.scatter(ch1_plot, ch2_plot, alpha=0.3, s=1, c='blue')
        ax5.set_xlabel('Channel 1 Intensity')
        ax5.set_ylabel('Channel 2 Intensity')
        ax5.set_title(f"Scatter Plot\nr={self.results['pearson_r']:.3f}",
                     fontsize=12, fontweight='bold')
        ax5.grid(True, alpha=0.3)

        ax6 = plt.subplot(3, 4, 6)
        h = ax6.hist2d(ch1_flat, ch2_flat, bins=50, cmap='viridis',
                       cmin=1, norm=LogNorm())
        plt.colorbar(h[3], ax=ax6, label='Count (log)')
        ax6.set_xlabel('Channel 1 Intensity')
        ax6.set_ylabel('Channel 2 Intensity')
        ax6.set_title('2D Histogram', fontsize=12, fontweight='bold')

        ax7 = plt.subplot(3, 4, 7)
        mid_row = self.ch1_raw.shape[0] // 2
        ax7.plot(self.ch1_raw[mid_row, :], 'r-', label='Ch1', linewidth=2)
        ax7.plot(self.ch2_raw[mid_row, :], 'g-', label='Ch2', linewidth=2)
        ax7.set_xlabel('X Position')
        ax7.set_ylabel('Intensity')
        ax7.set_title('Intensity Profile (mid-row)', fontsize=12, fontweight='bold')
        ax7.legend()
        ax7.grid(True, alpha=0.3)

        ax8 = plt.subplot(3, 4, 8)
        ax8.hist(self.ch1_raw[self.mask], bins=50, alpha=0.6,
                color='red', label='Ch1', density=True)
        ax8.hist(self.ch2_raw[self.mask], bins=50, alpha=0.6,
                color='green', label='Ch2', density=True)
        ax8.set_xlabel('Intensity')
        ax8.set_ylabel('Density')
        ax8.set_title('Intensity Distributions', fontsize=12, fontweight='bold')
        ax8.legend()
        ax8.grid(True, alpha=0.3)

        ax9 = plt.subplot(3, 4, 9)
        ax9.axis('off')
        metrics_text = f"""
        COLOCALIZATION METRICS
        {'='*35}

        Pearson's r: {self.results['pearson_r']:.4f}
        p-value: {self.results['pearson_p']:.4e}

        Manders (original):
          M1 = {self.results['manders_M1_original']:.4f}
          M2 = {self.results['manders_M2_original']:.4f}
        
        Manders (thresholded):
          tM1 = {self.results['manders_M1']:.4f}
          tM2 = {self.results['manders_M2']:.4f}

        Overlap coef: {self.results['overlap_coef']:.4f}

        Li's ICQ: {self.results['li_icq']:.4f}
        
        Cyto A = {self.results['cyto_a']:.4f}
        Cyto B = {self.results['cyto_b']:.4f}

        Colocalized pixels: {coloc_mask.sum()}
        Total pixels: {self.mask.sum()}
        Coloc fraction: {coloc_mask.sum()/self.mask.sum():.3f}
        """
        ax9.text(0.1, 0.5, metrics_text, fontsize=10,
                family='monospace', verticalalignment='center')

        ax10 = plt.subplot(3, 4, 10)
        ax10.imshow(self.ch1_thresh, cmap='Reds')
        ax10.set_title('Ch1 Thresholded', fontsize=12, fontweight='bold')
        ax10.axis('off')

        ax11 = plt.subplot(3, 4, 11)
        ax11.imshow(self.ch2_thresh, cmap='Greens')
        ax11.set_title('Ch2 Thresholded', fontsize=12, fontweight='bold')
        ax11.axis('off')

        ax12 = plt.subplot(3, 4, 12)
        coloc_intensity = np.zeros_like(self.ch1_raw)
        coloc_intensity[coloc_mask] = (self.ch1_raw[coloc_mask] +
                                       self.ch2_raw[coloc_mask]) / 2
        im = ax12.imshow(coloc_intensity, cmap='hot')
        #plt.colorbar(im, ax=ax12, label='Average Intensity')
        ax12.set_title('Colocalization Intensity', fontsize=12, fontweight='bold')
        ax12.axis('off')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"\nCombined figure saved to: {save_path}")

        plt.show()

        if save_individual and save_path:
            self._save_individual_plots(output_dir, ch1_flat, ch2_flat,
                                       ch1_plot, ch2_plot, merged,
                                       coloc_mask, coloc_intensity)

def _save_individual_plots(self, output_dir, ch1_flat, ch2_flat,
                           ch1_plot, ch2_plot, merged,
                           coloc_mask, coloc_intensity):

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(self.ch1_raw, cmap='Reds')
    ax.set_title('Channel 1', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / '01_channel1.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(self.ch2_raw, cmap='Greens')
    ax.set_title('Channel 2', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / '02_channel2.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(merged)
    ax.set_title('Merged Channels', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / '03_merged.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(coloc_mask, cmap='hot')
    ax.set_title('Colocalized Pixels', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / '04_colocalized_pixels.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(ch1_plot, ch2_plot, alpha=0.3, s=2, c='blue')
    ax.set_xlabel('Channel 1 Intensity', fontsize=12)
    ax.set_ylabel('Channel 2 Intensity', fontsize=12)
    ax.set_title(f"Scatter Plot (r={self.results['pearson_r']:.3f})",
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / '05_scatter_plot.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    h = ax.hist2d(ch1_flat, ch2_flat, bins=50, cmap='viridis',
                 cmin=1, norm=LogNorm())
    plt.colorbar(h[3], ax=ax, label='Count (log)')
    ax.set_xlabel('Channel 1 Intensity', fontsize=12)
    ax.set_ylabel('Channel 2 Intensity', fontsize=12)
    ax.set_title('2D Histogram', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / '06_2d_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    mid_row = self.ch1_raw.shape[0] // 2
    ax.plot(self.ch1_raw[mid_row, :], 'r-', label='Ch1', linewidth=1)
    ax.plot(self.ch2_raw[mid_row, :], 'g-', label='Ch2', linewidth=1)
    ax.set_xlabel('X Position', fontsize=12)
    ax.set_ylabel('Intensity', fontsize=12)
    ax.set_title('Intensity Profile (mid-row)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / '07_intensity_profile.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(self.ch1_raw[self.mask], bins=50, alpha=0.6,
           color='red', label='Ch1', density=True)
    ax.hist(self.ch2_raw[self.mask], bins=50, alpha=0.6,
           color='green', label='Ch2', density=True)
    ax.set_xlabel('Intensity', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Intensity Distributions', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / '08_intensity_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(self.ch1_thresh, cmap='Reds')
    ax.set_title('Channel 1 Thresholded', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / '09_ch1_thresholded.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(self.ch2_thresh, cmap='Greens')
    ax.set_title('Channel 2 Thresholded', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / '10_ch2_thresholded.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(coloc_intensity, cmap='hot')
  #  plt.colorbar(im, ax=ax, label='Average Intensity')
    ax.set_title('Colocalization Intensity', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / '11_colocalization_intensity.png', dpi=300, bbox_inches='tight')
    plt.close()

    metrics_text = f"""COLOCALIZATION ANALYSIS RESULTS
    {'='*50}
    
    Pearson's Correlation Coefficient:
      r = {self.results['pearson_r']:.4f}
      p-value = {self.results['pearson_p']:.4e}
      
    Manders (original):
    M1 = {self.results['manders_M1_original']:.4f}
    M2 = {self.results['manders_M2_original']:.4f}
    
    Manders' Coefficients:
      M1 (Ch1 overlapping with Ch2) = {self.results['manders_M1']:.4f}
      M2 (Ch2 overlapping with Ch1) = {self.results['manders_M2']:.4f}
    
    Overlap Coefficient = {self.results['overlap_coef']:.4f}
    
    Li's ICQ = {self.results['li_icq']:.4f}
      (Range: -0.5=segregated, 0=random, +0.5=colocalized)
    
    Cytofluorogram parameters:
    a (slope) = {self.results['cyto_a']:.4f}
    b (intercept) = {self.results['cyto_b']:.4f}
    
    Spatial Statistics:
      Colocalized pixels: {coloc_mask.sum()}
      Total pixels analyzed: {self.mask.sum()}
      Colocalization fraction: {coloc_mask.sum()/self.mask.sum():.4f}
    """
    with open(output_dir / 'metrics_summary.txt', 'w') as f:
        f.write(metrics_text)

    print('COLOCALIZATION ANALYSIS DONE')
