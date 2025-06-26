def predict_stability(ground_energy, vibrational_freqs):
    """
    Predicts if a molecule is stable.
    Args:
        ground_energy (float): Energy in Hartree
        vibrational_freqs (list): Frequencies in cm⁻¹
    Returns:
        str: "Stable" or "Unstable"
    """
    if any(freq < 0 for freq in vibrational_freqs):
        return "Unstable (imaginary frequencies)"
    return "Thermodynamically Stable"