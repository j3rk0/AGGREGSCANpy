import numpy as np

def calculate_a3v(sequence):
    """Calculate the a3v value for each amino acid based on propensities."""
    aa_propensities = {
    'I': 1.822, 'F': 1.754, 'V': 1.594, 'L': 1.38, 'Y': 1.159, 'W': 1.037, 'M': 0.91, 
    'C': 0.604, 'A': -0.036, 'T': -0.159, 'S': -0.294, 'P': -0.334, 'G': -0.535, 
    'K': -0.931, 'H': -1.033, 'Q': -1.231, 'R': -1.24, 'N': -1.302, 'E': -1.412, 'D': -1.836 }
    return np.array([aa_propensities[aa] for aa in sequence])

def calculate_a4v(a3v_sequence):

    # assign window size according to sequence length
    if a3v_sequence.shape[0] <= 75:
        window_size = 5
    elif a3v_sequence.shape[0] <= 175:
        window_size = 7
    elif a3v_sequence.shape[0] <= 300:
        window_size = 9
    else:
        window_size = 11

    # add virtual residues to the sequence
    virtual_a3v = np.concatenate(([-1.625],a3v_sequence,[-1.085]))

    """Calculate the a4v using a sliding window."""
    ret= np.convolve(virtual_a3v, np.ones(window_size)/window_size, mode='valid')
    start = [ret[0]] * (window_size//2 - 1)
    end = [ret[-1]] * (window_size//2 - 1)
    return np.concatenate((start,ret,end))


def identify_hot_spots(a4v, sequence):
    """Identify Hot Spots in the sequence where a4v > HST and no proline is present."""
    hst = -0.02 # threshold precomputed from swissprot
    hot_spots = []
    start = None
    
    for i, (value, aa) in enumerate(zip(a4v, sequence)):
        if value > hst and aa != 'P':
            if start is None:
                start = i
        else:
            if start is not None and i - start >= 5:  # Minimum 5 residues
                hot_spots.append((start, i-1))
            start = None
    return hot_spots



# Example of aggregating the profile for a sequence
def aggregate_profile(sequence):
    """Main function to compute the aggregation profile."""
    hst = -0.02 # threshold precomputed from swissprot
    a3v = calculate_a3v(sequence)
    a4v = calculate_a4v(a3v)
    hot_spots = identify_hot_spots(a4v, sequence)
    

    centered_a4v = a4v - hst

    # Calculate areas and metrics using trapezoidal integration
    TA = np.trapz(centered_a4v)

    a4v_above_tr = np.maximum(0,centered_a4v)

    AAT = np.trapz(a4v_above_tr)  
    
    # Hot Spot Area Calculation
    HSAs = []
    for start, end in hot_spots:
        HSAs.append( np.trapz(centered_a4v[start:end]))
    THSA = sum(HSAs)


    A4VHS = []
    for start, end in hot_spots:
        A4VHS.append( np.mean(a4v[start:end]))
    
    nhs = len(hot_spots) 


    return {
        'HST':hst,
        'a3vSA':np.mean(a3v),   
        'nHS':  nhs,
        'NnHS': nhs / len(sequence) * 100, 
        "AAT": AAT, # <- error
        "THSA": THSA,
        "TA": TA,
        'AAATr':AAT/len(sequence),
        'THSAr':THSA/len(sequence),
        'Na4vSS' : a4v.sum() / len(sequence) * 100,
        "Hot Spots": hot_spots,
        'HSAs':HSAs,
        'A4VHS':A4VHS,
        'a3v':a3v,
        'a4v':a4v
    }
