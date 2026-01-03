from skimage import io

def load_multichannel_tiff(filepath, ch1_idx=0, ch2_idx=1):
    img = io.imread(filepath)

    if img.ndim == 2:
        raise ValueError("Insert TIFF img.")

    elif img.ndim == 3:
        if img.shape[0] <= img.shape[2]:
            ch1 = img[ch1_idx]
            ch2 = img[ch2_idx]
        else:
            ch1 = img[:, :, ch1_idx]
            ch2 = img[:, :, ch2_idx]

    elif img.ndim == 4:
        ch1 = img[0, ch1_idx]
        ch2 = img[0, ch2_idx]

    else:
        raise ValueError("Invalid img dims")

    return ch1, ch2
