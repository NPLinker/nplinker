import math

# code to normalise peaks for spectral matching ("rosetta stone" stuff)
def sqrt_normalise(peaks):
    temp = []
    total = 0.0
    for mz, intensity in peaks:
        temp.append((mz, math.sqrt(intensity)))
        total += intensity
    norm_facc = math.sqrt(total)
    normalised_peaks = []
    for mz, intensity in temp:
        normalised_peaks.append((mz, intensity / norm_facc))
    return normalised_peaks
