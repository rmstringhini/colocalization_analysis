import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from skimage import io, filters, morphology
from skimage.util import img_as_float
from matplotlib.colors import LogNorm
from plotting import (
    plot_analysis,
    plot_cytofluorogram,
    _save_individual_plots
)    


class ColocalizationAnalysis:
    def __init__(self, ch1, ch2, threshold_method='Otsu'):
        self.ch1_raw = img_as_float(ch1)
        self.ch2_raw = img_as_float(ch2)
        self.threshold_method = threshold_method
        self.ch1_thresh = None
        self.ch2_thresh = None
        self.mask = None
        self.results = {}
        
    def set_fiji_mask(self, include_zero_zero=False):
        if include_zero_zero:
            self.mask = np.ones_like(self.ch1_raw, dtype=bool)
        else:
            self.mask = (self.ch1_raw > 0) | (self.ch2_raw > 0)
    
    def apply_threshold(self, ch1_manual=None, ch2_manual=None):
        if self.threshold_method == 'Otsu':
            t1 = filters.threshold_otsu(self.ch1_raw)
            t2 = filters.threshold_otsu(self.ch2_raw)
            self.ch1_thresh = self.ch1_raw > t1
            self.ch2_thresh = self.ch2_raw > t2
            print("Otsu threshold applied")
            # print(f"Otsu threshold Ch1: {t1}")
            # print(f"Otsu threshold Ch2: {t2}")
            
        elif self.threshold_method == 'manual':
            if ch1_manual is None or ch2_manual is None:
                raise ValueError("Manual thresholding required")
            self.ch1_thresh = self.ch1_raw > ch1_manual
            self.ch2_thresh = self.ch2_raw > ch2_manual
            print("Manual threshold applied")
        
        else:
            self.ch1_thresh = self.ch1_raw > 0
            self.ch2_thresh = self.ch2_raw > 0
            
        self.mask = self.ch1_thresh | self.ch2_thresh
    
    def pearson_correlation(self):
        ch1_masked = self.ch1_raw[self.mask]
        ch2_masked = self.ch2_raw[self.mask]
    
        r, p_value = stats.pearsonr(ch1_masked, ch2_masked)
        self.results['pearson_r'] = r
        self.results['pearson_p'] = p_value
        
        return r, p_value
        
    def manders_original(self):
        ch1 = self.ch1_raw
        ch2 = self.ch2_raw
    
        M1 = ch1[ch2 > 0].sum() / ch1.sum()
        M2 = ch2[ch1 > 0].sum() / ch2.sum()
    
        self.results['manders_M1_original'] = M1
        self.results['manders_M2_original'] = M2
    
        return M1, M2
        
    def manders_coefficients(self):
        ch1_masked = self.ch1_raw[self.mask]
        ch2_masked = self.ch2_raw[self.mask]
        
        coloc_mask = self.ch1_thresh & self.ch2_thresh
        
        ch1_coloc = self.ch1_raw[coloc_mask].sum()
        ch2_coloc = self.ch2_raw[coloc_mask].sum()
        
        ch1_total = self.ch1_raw[self.ch1_thresh].sum()
        ch2_total = self.ch2_raw[self.ch2_thresh].sum()
        
        M1 = ch1_coloc / ch1_total if ch1_total > 0 else 0
        M2 = ch2_coloc / ch2_total if ch2_total > 0 else 0
        
        self.results['manders_M1'] = M1
        self.results['manders_M2'] = M2
        
        return M1, M2
    
    
    def overlap_coefficient(self):
        ch1_masked = self.ch1_raw[self.mask]
        ch2_masked = self.ch2_raw[self.mask]
        
        numerator = np.sum(ch1_masked * ch2_masked)
        denominator = np.sqrt(np.sum(ch1_masked**2) * np.sum(ch2_masked**2))
        
        r = numerator / denominator if denominator > 0 else 0
        self.results['overlap_coef'] = r
        
        return r
    
    def costes_threshold(self):
        from scipy.optimize import minimize_scalar
        from scipy.stats import pearsonr
        
        def negative_correlation(t):
            mask = (self.ch1_raw > t) | (self.ch2_raw > t)
            if mask.sum() < 10:
                return 1.0
            ch1_m = self.ch1_raw[mask]
            ch2_m = self.ch2_raw[mask]
            try:
                r, _ = stats.pearsonr(ch1_m, ch2_m)
                return abs(r)
            except:
                return 1.0
        
        res1 = minimize_scalar(negative_correlation, bounds=(0,1), method='bounded')
        return res1.x
    
    def li_icq(self):
        ch1_masked = self.ch1_raw[self.mask]
        ch2_masked = self.ch2_raw[self.mask]
        
        mean_ch1 = np.mean(ch1_masked)
        mean_ch2 = np.mean(ch2_masked)
        
        pdm = (ch1_masked - mean_ch1) * (ch2_masked - mean_ch2)
        
        num_positive = np.sum(pdm > 0)
        total = len(pdm)
        
        icq = (num_positive / total) - 0.5
        self.results['li_icq'] = icq
        return icq
    
    def cytofluorogram_regression(self):    
        x = self.ch1_raw[self.mask].astype(np.float64)
        y = self.ch2_raw[self.mask].astype(np.float64)
    
        # Means
        x_mean = x.mean()
        y_mean = y.mean()
    
        # Variance and covariance
        cov_xy = np.mean((x - x_mean) * (y - y_mean))
        var_x  = np.mean((x - x_mean) ** 2)
    
        a = cov_xy / var_x if var_x != 0 else 0
        b = y_mean - a * x_mean
        b_aux = b
        b = b * 255
        
        self.results['cyto_a'] = a
        self.results['cyto_b'] = b
    
        return a, b, b_aux
    
    def compute_cyto_regression(self):
        ch1 = self.ch1_raw.flatten()
        ch2 = self.ch2_raw.flatten()
    
        valid = (ch1 + ch2) > 0
        ch1 = ch1[valid]
        ch2 = ch2[valid]

        mean_x = np.mean(ch1)
        mean_y = np.mean(ch2)
    
        cov = np.mean((ch1 - mean_x) * (ch2 - mean_y))
        var = np.mean((ch1 - mean_x) ** 2)
    
        a = cov / var
        b = mean_y - a * mean_x
    
        return a, b

    

    def run_analysis(self, ch1_manual=None, ch2_manual=None):
        print("============================================")
        print('COLOCALIZATION ANALYSIS')
        
        self.apply_threshold(ch1_manual, ch2_manual)
        
        self.set_fiji_mask(include_zero_zero=True)
        
        print('METRICS')
        r, p = self.pearson_correlation()
        print('Pearson R:', r, 'Pearson P:', p)
        
        M1_orig, M2_orig = self.manders_original()
        print('Manders (original)')
        print('M1:', M1_orig, 'M2:', M2_orig)
        
        M1, M2 = self.manders_coefficients()
        print('Manders (thresholded)')
        print('M1:', M1, 'M2:', M2)
        
        oc = self.overlap_coefficient()
        print('Overlap Coefficient:', oc)
        
        icq = self.li_icq()
        print('ICQ:', icq)
        
        a, b, b_aux = self.cytofluorogram_regression()
        print('Cytofluorogram parameters:', 'a: ',a , 'b: ',b)
        
        
        # a_cyto_plot, b_cyto_plot = self.compute_cyto_regression()
        a_cyto_plot, b_cyto_plot, b_aux = self.cytofluorogram_regression()

        print(f"Cytofluorogram parameters:")
        print(f"a (slope): {a_cyto_plot}")
        print(f"b (intercept): {b_cyto_plot}")
        print("============================================")
        
        self.plot_cytofluorogram(
            a=a_cyto_plot,
            b=b_aux,
            ch1_thresh=np.mean(self.ch1_raw[self.ch1_thresh]),
            ch2_thresh=np.mean(self.ch2_raw[self.ch2_thresh]),
            save_path="D:/colocalization_test_img/colocalization_analysis/data_and_results/cytofluorogram.png"
            )
        
        return self.results
    
ColocalizationAnalysis.plot_analysis = plot_analysis
ColocalizationAnalysis.plot_cytofluorogram = plot_cytofluorogram
ColocalizationAnalysis._save_individual_plots = _save_individual_plots