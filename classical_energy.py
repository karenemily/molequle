from pyscf import gto, scf, hessian

def calculate_ground_state_energy(molecule, basis="sto-3g"):
    """
    Compute ground-state energy using Hartree-Fock (PySCF).
    Args:
        molecule (str): e.g., "H 0 0 0; F 0 0 1.1"
        basis (str): Basis set (e.g., "sto-3g", "cc-pvdz")
    Returns:
        float: Ground-state energy (in Hartree)
    """
    mol = gto.M(atom=molecule, basis=basis)
    hf = scf.RHF(mol).run()
    return hf.e_tot

def compute_vibrational_frequencies(molecule):
    """
    Check stability via vibrational frequencies.
    Returns:
        list: Frequencies (cm⁻¹). If any are negative, molecule is unstable.
    """
    mol = gto.M(atom=molecule, basis="sto-3g")
    hf = scf.RHF(mol).run()
    hess = hessian.RHF(hf).run()
    return hess.vib_freq()