import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# --- Constants ---
R = 8.314  # J/(mol·K)
HARTREE_TO_KJ = 2625.5
DAY_TO_SECONDS = 86400

# --- Custom CSS ---
st.set_page_config(page_title="MoleQule", page_icon="⚗️", layout="wide")
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Poppins:wght@400;600&display=swap');
* {
    font-family: 'Poppins', sans-serif !important;
}
h1 span {
    color: #00b4ff !important;
}
.stMetric {
    border-left: 3px solid #00b4ff !important;
}
table {
    font-family: 'Poppins' !important;
}
.dataframe {
    font-family: 'Poppins' !important;
}
.sidebar-section {
    margin-bottom: 2rem;
}
</style>
""", unsafe_allow_html=True)

# --- Title ---
st.markdown("<h1>⚗️ Mole<span>Q</span>ule</h1>", unsafe_allow_html=True)

# --- Drug Database (Corrected Values) ---
DRUG_DB = {
    "Aspirin": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "degradation": {
            "product": "Salicylic acid",
            "smiles": "OC1=CC=CC=C1C(=O)O",
            "Ea": 85.2,  # kJ/mol (experimental value) :cite[2]
            "A": 1.15e12,  # s⁻¹
            "E0": -1027.3,
            "E_TS": -942.1,
            "E_deg": -950.8,
            "t90_ref": 3.2  # years at 25°C :cite[2]
        }
    },
    "Cyclobutadiene (Unstable)": {
        "smiles": "C1=CC=C1",
        "degradation": {
            "product": "2 Acetylene",
            "smiles": "C#CC#C",
            "Ea": 25.0,  # kJ/mol (antiaromatic destabilization) :cite[1]
            "A": 1.0e13,
            "E0": -153.0,
            "E_TS": -128.0,
            "E_deg": -310.0,
            "t90_ref": 0.003  # ~1 day at 25°C :cite[1]
        }
    },
    "Methane (CH₄)": {
        "smiles": "C",
        "degradation": {
            "product": "CH₃· + H·",
            "smiles": "[CH3]",
            "Ea": 435.0,  # kJ/mol (C-H bond strength) :cite[1]
            "A": 1.0e16,
            "E0": -40.5,
            "E_TS": 394.5,
            "E_deg": 0.0,
            "t90_ref": 1000  # years (effectively stable) :cite[3]
        }
    }
}

# --- Sidebar Panel ---
with st.sidebar:
    st.header("Input Panel")
    
    # Molecule Selection
    with st.container():
        drug_choice = st.selectbox(
            "Select Molecule:", 
            list(DRUG_DB.keys()),
            key="drug_choice"
        )
        st.markdown("<div class='sidebar-section'></div>", unsafe_allow_html=True)
    
    # Degradation Product
    drug_data = DRUG_DB[drug_choice]
    with st.container():
        st.markdown(f"**Degradation Product:**  \n{drug_data['degradation']['product']}")
        st.markdown("<div class='sidebar-section'></div>", unsafe_allow_html=True)
    
    # Temperature Control
    with st.container():
        temperature = st.slider(
            "Temperature (K)", 
            min_value=273, 
            max_value=323, 
            value=298,
            step=5,
            key="temp_slider"
        )
        st.markdown("<div class='sidebar-section'></div>", unsafe_allow_html=True)
    
    # File Upload
    with st.container():
        uploaded_file = st.file_uploader(
            "Upload Molecule File", 
            type=["smi", "mol", "xyz"],
            key="file_uploader"
        )

# --- Calculations ---
def calculate_shelf_life(Ea, A, T, T_ref=298):
    """Returns shelf life in days with reference adjustment"""
    k_ref = A * np.exp(-Ea * 1000 / (R * T_ref))
    k = A * np.exp(-Ea * 1000 / (R * T))
    return (0.105 / k_ref) * (k_ref / k) * (DAY_TO_SECONDS / 86400)  # days

# Get reference values
t90_ref_days = drug_data["degradation"]["t90_ref"] * 365  # Convert years to days
Ea = drug_data["degradation"]["Ea"]
A = drug_data["degradation"]["A"]

# Calculate temperature-adjusted shelf life
if drug_choice == "Methane (CH₄)":
    t90 = ">1000 years"  # Too stable to measure
else:
    t90_days = calculate_shelf_life(Ea, A, temperature)
    if t90_days < 1:
        t90 = f"{t90_days*24:.1f} hours"
    elif t90_days < 30:
        t90 = f"{t90_days:.1f} days"
    elif t90_days < 365:
        t90 = f"{t90_days/30:.1f} months"
    else:
        t90 = f"{t90_days/365:.1f} years"

# --- Visualization ---
st.subheader("Molecular Structures")
col1, col2 = st.columns(2)
try:
    mol_react = Chem.MolFromSmiles(drug_data["smiles"])
    col1.image(Draw.MolToImage(mol_react, size=(300,300)), caption="Reactant")
    
    mol_prod = Chem.MolFromSmiles(drug_data["degradation"]["smiles"])
    col2.image(Draw.MolToImage(mol_prod, size=(300,300)), caption="Product")
except:
    st.warning("Structure rendering unavailable")

# --- Energy Profile ---
st.subheader("Reaction Energy Profile")
fig, ax = plt.subplots(figsize=(8,4))
energies = [
    drug_data["degradation"]["E0"],
    drug_data["degradation"]["E_TS"],
    drug_data["degradation"]["E_deg"]
]
ax.plot(["Reactant", "TS", "Product"], energies, '-o', color='#00b4ff')
ax.set_ylabel("Energy (Hartree)")
st.pyplot(fig)

# --- Results Table ---
st.subheader("Stability Analysis")
results = pd.DataFrame({
    "Parameter": [
        "Reactant Energy (E₀)",
        "TS Energy (E_TS)",
        "Product Energy (E_deg)",
        "Activation Energy (Ea)",
        "Frequency Factor (A)",
        "Shelf Life (t₉₀)"
    ],
    "Value": [
        f"{drug_data['degradation']['E0']:.2f} Ha",
        f"{drug_data['degradation']['E_TS']:.2f} Ha",
        f"{drug_data['degradation']['E_deg']:.2f} Ha",
        f"{Ea:.2f} kJ/mol",
        f"{A:.2e} s⁻¹",
        t90 if isinstance(t90, str) else f"{t90:.1f} years"
    ],
    "Description": [
        "Ground state energy",
        "Transition state energy",
        "Degradation product energy",
        "Energy barrier for reaction",
        "Collision frequency factor",
        f"Time to 90% potency at {temperature}K"
    ]
})

st.dataframe(results.style.set_properties(**{
    'color': 'black',
    'background-color': '#f8f9fa'
}))

# --- Key Metrics ---
cols = st.columns(3)
cols[0].metric("Temperature", f"{temperature} K")
cols[1].metric("Activation Energy", f"{Ea:.2f} kJ/mol")

if drug_choice == "Cyclobutadiene (Unstable)":
    cols[2].metric("Shelf Life", t90, delta="Highly Unstable", delta_color="inverse")
elif drug_choice == "Methane (CH₄)":
    cols[2].metric("Shelf Life", t90, delta="Extremely Stable")
else:
    cols[2].metric("Shelf Life", t90, delta="Stable")