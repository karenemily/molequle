from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_algorithms import VQE
from qiskit_aer.primitives import Estimator
from qiskit.circuit.library import EfficientSU2

def calculate_quantum_energy(molecule, basis="sto-3g"):
    """
    Compute ground-state energy using VQE (Qiskit).
    Args:
        molecule (str): e.g., "H 0 0 0; H 0 0 0.74"
    Returns:
        float: Ground-state energy (in Hartree)
    """
    driver = PySCFDriver(atom=molecule, basis=basis)
    problem = driver.run()
    mapper = JordanWignerMapper()
    ansatz = EfficientSU2(problem.num_spatial_orbitals)
    vqe = VQE(Estimator(), ansatz)
    solver = GroundStateEigensolver(mapper, vqe)
    result = solver.solve(problem)
    return result.total_energies[0]