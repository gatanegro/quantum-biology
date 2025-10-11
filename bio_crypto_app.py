import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import time
from scipy import stats
"""
LOGOS THEORY
Author: Martin Doina 
"""

# Constants for Bridge Formula
LZ = 1.23488369648610768936904233243207400555883244095073 
QDF = 0.8097928609310675
me = 9.10938356e-31  # kg
pi = np.pi
c = 299792458  # m/s

# Cryptochrome domain references
crypto_mass_ref = 2.5e-25  # kg (flavin molecule scale)
crypto_energy_ref = 2.0e-19  # J (magnetic field interactions)  
crypto_radius_ref = 2.0e-9  # m (2 nm - protein domain spacing)

def F_HQS_tanh(n):
    """Harmonic Quantum Stability - tanh resonance"""
    HQS3 = 0.235543306845346130248413363583896203088885727325439
    return 1.0 + 0.01 * np.tanh(n * HQS3)

def compute_spatial_radius(n, a0_prime=None):
    """Spatial formula for cryptochrome systems"""
    if a0_prime is None:
        a0_prime = crypto_radius_ref
    
    resonance_factor = F_HQS_tanh(n)
    radius = a0_prime * (LZ ** n) * resonance_factor
    return radius, resonance_factor

def compute_n_mass_crypto(mass):
    """Bridge Formula with cryptochrome mass reference"""
    return pi * (np.log10(mass / crypto_mass_ref) - np.log10(QDF)) / np.log10(LZ)

def octave_position(n):
    """Map n-value to circular octave position"""
    reduced_n = int(round(n)) % 9
    angle = 2 * pi * reduced_n / 9
    return np.cos(angle), np.sin(angle), reduced_n, angle

def get_mass_by_resname_crypto(resname):
    """Cryptochrome and magnetoreception molecule database"""
    crypto_masses = {
        # Flavin molecules (key for magnetoreception)
        'FAD': 785.55, 'FMN': 456.34, 'FMN': 456.34,
        'RBF': 785.55,  # Riboflavin
        # Tryptophan radicals (quantum coherence)
        'TRP': 204.23, 'W': 204.23,
        # Electron transfer molecules
        'HEM': 616.49, 'HEC': 616.49, 'HEA': 616.49,
        'CYS': 121.16, 'TYR': 181.19, 'HIS': 155.16,
        # Metal centers
        'ZN': 65.38, 'FE': 55.85, 'CU': 63.55, 'MG': 24.31,
        # Cryptochrome specific
        'CR1': 500.0, 'CR2': 500.0, 'CRY': 500.0,
        # Other redox-active residues
        'MET': 149.21, 'ARG': 174.20, 'LYS': 146.19,
    }
    amu_to_kg = 1.66053906660e-27
    return crypto_masses.get(resname, 300.0) * amu_to_kg  # Default for small molecules

def extract_crypto_data_pdb(pdb_file):
    """Extract magnetoreception-related molecules from PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    crypto_data = []
    
    crypto_molecules = ['FAD', 'FMN', 'RBF', 'TRP', 'W', 'HEM', 'HEC', 'HEA', 
                       'CYS', 'TYR', 'HIS', 'ZN', 'FE', 'CU', 'MG', 'CR1', 'CR2', 'CRY']
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in crypto_molecules:
                    mass = get_mass_by_resname_crypto(resname)
                    crypto_data.append({
                        'resname': resname,
                        'mass': mass,
                        'source': pdb_file,
                        'type': 'flavin' if resname in ['FAD','FMN','RBF'] else 
                               'tryptophan' if resname in ['TRP','W'] else
                               'metal' if resname in ['ZN','FE','CU','MG'] else
                               'redox' if resname in ['CYS','TYR','HIS'] else
                               'other'
                    })
    
    return crypto_data

def extract_crypto_data_cif(cif_file):
    """Extract cryptochrome molecules from CIF file"""
    d = MMCIF2Dict(cif_file)
    crypto_data = []
    
    resnames = d.get('_atom_site.label_comp_id', [])
    if isinstance(resnames, str):
        resnames = [resnames]
    
    crypto_molecules = ['FAD', 'FMN', 'RBF', 'TRP', 'W', 'HEM', 'HEC', 'HEA', 
                       'CYS', 'TYR', 'HIS', 'ZN', 'FE', 'CU', 'MG']
    
    residue_dict = {}
    for resname in resnames:
        resname_clean = resname.strip()
        if resname_clean in crypto_molecules and resname_clean not in residue_dict:
            mass = get_mass_by_resname_crypto(resname_clean)
            residue_dict[resname_clean] = {
                'resname': resname_clean,
                'mass': mass,
                'source': cif_file,
                'type': 'flavin' if resname_clean in ['FAD','FMN','RBF'] else 
                       'tryptophan' if resname_clean in ['TRP','W'] else
                       'metal' if resname_clean in ['ZN','FE','CU','MG'] else
                       'redox'
            }
    
    crypto_data = list(residue_dict.values())
    return crypto_data

def calculate_quantum_metrics_crypto(crypto_data):
    """Calculate quantum coherence metrics for cryptochrome"""
    if not crypto_data:
        return {}
    
    masses = [p['mass'] for p in crypto_data]
    n_values = [compute_n_mass_crypto(mass) for mass in masses]
    octaves = [octave_position(n)[2] for n in n_values]
    
    mass_stats = {
        'mean': np.mean(masses),
        'std': np.std(masses),
        'cv': (np.std(masses) / np.mean(masses)) * 100,
        'range': (min(masses), max(masses))
    }
    
    n_stats = {
        'mean': np.mean(n_values),
        'std': np.std(n_values),
        'range': (min(n_values), max(n_values))
    }
    
    octave_counts = {}
    for octave in octaves:
        octave_counts[octave] = octave_counts.get(octave, 0) + 1
    
    superposition_octaves = [octave for octave, count in octave_counts.items() if count > 1]
    superposition_ratio = len([p for p in crypto_data if octave_counts[octave_position(compute_n_mass_crypto(p['mass']))[2]] > 1]) / len(crypto_data)
    
    # Group by molecule type
    type_counts = {}
    for item in crypto_data:
        mol_type = item['type']
        type_counts[mol_type] = type_counts.get(mol_type, 0) + 1
    
    return {
        'molecule_count': len(crypto_data),
        'type_distribution': type_counts,
        'mass_stats': mass_stats,
        'n_stats': n_stats,
        'octave_distribution': octave_counts,
        'superposition_octaves': superposition_octaves,
        'superposition_ratio': superposition_ratio,
        'mass_tuning': mass_stats['cv'],
        'flavin_present': 'flavin' in type_counts
    }

# Global data storage
all_crypto_data = []

def load_crypto_pdb():
    file_path = filedialog.askopenfilename(filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
    if not file_path:
        return
    crypto_molecules = extract_crypto_data_pdb(file_path)
    if crypto_molecules:
        all_crypto_data.extend(crypto_molecules)
        update_crypto_plot()
        update_crypto_analysis(crypto_molecules, "PDB")
    else:
        messagebox.showwarning("No Molecules", "No cryptochrome-related molecules found.")

def load_crypto_cif():
    file_path = filedialog.askopenfilename(filetypes=[("CIF files", "*.cif *.mmcif"), ("All files", "*.*")])
    if not file_path:
        return
    crypto_molecules = extract_crypto_data_cif(file_path)
    if crypto_molecules:
        all_crypto_data.extend(crypto_molecules)
        update_crypto_plot()
        update_crypto_analysis(crypto_molecules, "CIF")
    else:
        messagebox.showwarning("No Molecules", "No cryptochrome-related molecules found.")

def update_crypto_analysis(crypto_data, file_type):
    analysis_text.delete(1.0, tk.END)
    
    try:
        metrics = calculate_quantum_metrics_crypto(crypto_data)
        
        text = "MAGNETORECEPTION QUANTUM ANALYSIS\n"
        text += "=" * 50 + "\n\n"
        
        text += f"STRUCTURE: {file_type}\n"
        text += f"Molecules analyzed: {metrics['molecule_count']}\n"
        text += f"Molecule types: {metrics['type_distribution']}\n\n"
        
        text += "QUANTUM STATE ANALYSIS:\n"
        text += f"  Recursive depth: {metrics['n_stats']['mean']:.1f} ± {metrics['n_stats']['std']:.1f} cycles\n"
        text += f"  Current phases: {dict(metrics['octave_distribution'])}\n"
        text += f"  Mass coherence: {metrics['mass_tuning']:.2f}%\n\n"
        
        text += "QUANTUM COHERENCE:\n"
        text += f"  Superposition ratio: {metrics['superposition_ratio']:.1%}\n"
        if metrics['superposition_octaves']:
            text += f"  Coherent phases: {metrics['superposition_octaves']}\n\n"
        else:
            text += "  No superposition detected\n\n"
        
        # Magnetic sensitivity prediction
        if metrics['flavin_present'] and metrics['superposition_ratio'] > 0.3:
            text += "MAGNETIC SENSITIVITY: HIGH\n"
            text += "  • Flavin radicals present ✓\n"
            text += "  • Quantum coherence detected ✓\n"
            text += "  • Likely magnetoreception capable\n\n"
        elif metrics['flavin_present']:
            text += "MAGNETIC SENSITIVITY: MEDIUM\n"
            text += "  • Flavin radicals present ✓\n"
            text += "  • Limited quantum coherence\n\n"
        else:
            text += "MAGNETIC SENSITIVITY: LOW\n"
            text += "  • No flavin radicals detected\n\n"
        
        # Spatial organization
        if crypto_data:
            n_values = [compute_n_mass_crypto(p['mass']) for p in crypto_data]
            avg_n = np.mean(n_values)
            spacing, resonance = compute_spatial_radius(avg_n)
            text += f"SPATIAL ORGANIZATION:\n"
            text += f"  Optimal spacing: {spacing*1e9:.1f} nm\n"
            text += f"  Resonance factor: {resonance:.3f}\n\n"
        
        n_mean = metrics['n_stats']['mean']
        completed_cycles = int(n_mean)
        current_phase = n_mean - completed_cycles
        
        text += "QUANTUM INTERPRETATION:\n"
        text += f"  Evolutionary depth: ~{completed_cycles} recursive cycles\n"
        text += f"  Current cycle: {current_phase:.1%} complete\n"
        text += f"  Active quantum phases: {len(metrics['octave_distribution'])}/9\n"
        if metrics['superposition_octaves']:
            text += f"  Phase {metrics['superposition_octaves'][0]} shows quantum coherence\n"
        
        analysis_text.insert(1.0, text)
        
    except Exception as e:
        analysis_text.insert(1.0, f"Error in analysis: {str(e)}")

def update_crypto_plot():
    ax.cla()
    
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.set_aspect('equal')
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    
    for radius in [0.25, 0.5, 0.75, 1.0]:
        circle = plt.Circle((0, 0), radius, color='lightgray', fill=False, linestyle='--', alpha=0.5)
        ax.add_artist(circle)
    
    main_circle = plt.Circle((0, 0), 1.0, color='black', fill=False, linewidth=1.5)
    ax.add_artist(main_circle)
    
    for i in range(9):
        angle = 2 * pi * i / 9
        x = np.cos(angle)
        y = np.sin(angle)
        ax.plot([0, x*1.2], [0, y*1.2], 'gray', linestyle='--', alpha=0.7, linewidth=0.8)
        ax.text(x*1.25, y*1.25, f'O{i}', fontsize=9, color='darkblue')
    
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    
    if all_crypto_data:
        color_map = {
            'flavin': 'yellow', 'tryptophan': 'purple', 
            'metal': 'red', 'redox': 'blue', 'other': 'green'
        }
        
        for molecule in all_crypto_data:
            n_val = compute_n_mass_crypto(molecule['mass'])
            x, y, octave, angle = octave_position(n_val)
            color = color_map.get(molecule['type'], 'black')
            
            ax.plot(x, y, 'o', color=color, markersize=8, alpha=0.8, markeredgecolor='black', markeredgewidth=1)
            ax.text(x + 0.05, y + 0.05, molecule['resname'], fontsize=7)
    
    if all_crypto_data:
        n_values = [compute_n_mass_crypto(p['mass']) for p in all_crypto_data]
        avg_cycles = np.mean(n_values)
        completed_cycles = int(avg_cycles)
        ax.set_title(f"Cryptochrome Analysis - {completed_cycles} Cycles")
    else:
        ax.set_title("Cryptochrome Quantum Analysis")
    
    canvas.draw()

def clear_crypto_analysis():
    global all_crypto_data
    all_crypto_data = []
    analysis_text.delete(1.0, tk.END)
    update_crypto_plot()

# Create interface
root = tk.Tk()
root.title("Cryptochrome Quantum Analyzer")
root.geometry("1400x800")

# Main container
main_frame = ttk.Frame(root, padding="15")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Control panel
control_frame = ttk.Frame(main_frame, padding="10")
control_frame.grid(row=0, column=0, sticky=(tk.N, tk.S), padx=(0, 10))

ttk.Button(control_frame, text="Load Crypto PDB", command=load_crypto_pdb, width=15).grid(row=0, column=0, pady=5)
ttk.Button(control_frame, text="Load Crypto CIF", command=load_crypto_cif, width=15).grid(row=1, column=0, pady=5)
ttk.Button(control_frame, text="Clear", command=clear_crypto_analysis, width=15).grid(row=2, column=0, pady=5)

# Analysis display
analysis_text = tk.Text(control_frame, width=35, height=25, wrap=tk.WORD, font=('Courier', 9), bg='white')
analysis_text.grid(row=3, column=0, pady=10)

# Plot area
plot_frame = ttk.Frame(main_frame)
plot_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))

# Create figure
fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
fig.patch.set_facecolor('white')

canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# Configure grid weights
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
main_frame.columnconfigure(1, weight=1)
main_frame.rowconfigure(0, weight=1)

# Initialize
update_crypto_plot()

root.mainloop()
