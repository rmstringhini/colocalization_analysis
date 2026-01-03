from analysis import ColocalizationAnalysis
from utils import load_multichannel_tiff

def analyze_from_file(filepath, ch1_idx=0, ch2_idx=1,
                      threshold_method='Otsu',
                      save_results=None):

    ch1, ch2 = load_multichannel_tiff(filepath, ch1_idx, ch2_idx)

    coloc = ColocalizationAnalysis(ch1, ch2, threshold_method)

    coloc.run_analysis()

    if save_results is None:
        save_results = filepath.replace('.tif', '_colocalization.png')

    coloc.plot_analysis(save_path=save_results)

    return coloc
