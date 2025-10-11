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

# DNA/RNA domain references
dna_mass_ref = 3.0e-25  # kg (nucleotide scale)
dna_energy_ref = 4.0e-19  # J (UV damage repair range)
dna_radius_ref = 0.34e-9  # m (0.34 nm - base pair spacing)

def F_HQS_tanh(n):
    """Harmonic Quantum Stability - tanh resonance"""
    HQS3 = 0.235543306845346130248413363583896203088885727325439
    return 1.0 + 0.01 * np.tanh(n * HQS3)

def compute_spatial_radius(n, a0_prime=None):
    """Spatial formula for DNA systems"""
    if a0_prime is None:
        a0_prime = dna_radius_ref
    
    resonance_factor = F_HQS_tanh(n)
    radius = a0_prime * (LZ ** n) * resonance_factor
    return radius, resonance_factor

def compute_n_mass_dna(mass):
    """Bridge Formula with DNA mass reference"""
    return pi * (np.log10(mass / dna_mass_ref) - np.log10(QDF)) / np.log10(LZ)

def octave_position(n):
    """Map n-value to circular octave position"""
    reduced_n = int(round(n)) % 9
    angle = 2 * pi * reduced_n / 9
    return np.cos(angle), np.sin(angle), reduced_n, angle

def get_mass_by_resname_dna(resname):
    """DNA/RNA molecule database"""
    dna_masses = {
        # DNA bases
        'DA': 313.21, 'DT': 304.20, 'DC': 289.18, 'DG': 329.21,
        'ADE': 313.21, 'THY': 304.20, 'CYT': 289.18, 'GUA': 329.21,
        'A': 313.21, 'T': 304.20, 'C': 289.18, 'G': 329.21,
        # RNA bases  
        'RA': 313.21, 'RU': 306.17, 'RC': 289.18, 'RG': 329.21,
        'URA': 306.17, 'U': 306.17,
        # Modified bases
        '5MC': 303.20, 'H2U': 302.18, 'PSU': 324.20,
        # Backbone components
        'DPR': 97.07,  # Deoxyribose approximate
        'RIB': 115.09,  # Ribose approximate
        'PO4': 94.97,   # Phosphate group
        # Ions and cofactors
        'MG': 24.31, 'CA': 40.08, 'ZN': 65.38, 'MN': 54.94,
        'K': 39.10, 'NA': 22.99,
        # DNA repair enzymes components
        'FAD': 785.55, 'FMN': 456.34,  # Photolyase cofactors
        'HIS': 155.16, 'CYS': 121.16, 'TYR': 181.19,  # Catalytic residues
        # Water and hydration
        'HOH': 18.02, 'DOD': 20.03,
    }
    amu_to_kg = 1.66053906660e-27
    
    # Handle nucleic acids specifically
    if resname in ['DA', 'DT', 'DC', 'DG', 'A', 'T', 'C', 'G']:
        return dna_masses.get(resname, 300.0) * amu_to_kg
    elif resname in ['RA', 'RU', 'RC', 'RG', 'U']:
        return dna_masses.get(resname, 300.0) * amu_to_kg
    else:
        return dna_masses.get(resname, 200.0) * amu_to_kg

def extract_dna_data_pdb(pdb_file):
    """Extract DNA/RNA molecules from PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    dna_data = []
    
    dna_molecules = ['DA', 'DT', 'DC', 'DG', 'ADE', 'THY', 'CYT', 'GUA', 'A', 'T', 'C', 'G',
                    'RA', 'RU', 'RC', 'RG', 'URA', 'U', '5MC', 'H2U', 'PSU',
                    'DPR', 'RIB', 'PO4', 'MG', 'CA', 'ZN', 'MN', 'K', 'NA',
                    'FAD', 'FMN', 'HIS', 'CYS', 'TYR', 'HOH', 'DOD']
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in dna_molecules:
                    mass = get_mass_by_resname_dna(resname)
                    dna_data.append({
                        'resname': resname,
                        'mass': mass,
                        'source': pdb_file,
                        'type': 'dna_base' if resname in ['DA','DT','DC','DG','A','T','C','G'] else
                               'rna_base' if resname in ['RA','RU','RC','RG','U'] else
                               'modified_base' if resname in ['5MC','H2U','PSU'] else
                               'backbone' if resname in ['DPR','RIB','PO4'] else
                               'ion' if resname in ['MG','CA','ZN','MN','K','NA'] else
                               'cofactor' if resname in ['FAD','FMN'] else
                               'catalytic' if resname in ['HIS','CYS','TYR'] else
                               'water'
                    })
    
    return dna_data

def extract_dna_data_cif(cif_file):
    """Extract DNA/RNA molecules from CIF file"""
    d = MMCIF2Dict(cif_file)
    dna_data = []
    
    resnames = d.get('_atom_site.label_comp_id', [])
    if isinstance(resnames, str):
        resnames = [resnames]
    
    dna_molecules = ['DA', 'DT', 'DC', 'DG', 'A', 'T', 'C', 'G',
                    'RA', 'RU', 'RC', 'RG', 'U', 'MG', 'CA', 'ZN']
    
    residue_dict = {}
    for resname in resnames:
        resname_clean = resname.strip()
        if resname_clean in dna_molecules and resname_clean not in residue_dict:
            mass = get_mass_by_resname_dna(resname_clean)
            residue_dict[resname_clean] = {
                'resname': resname_clean,
                'mass': mass,
                'source': cif_file,
                'type': 'dna_base' if resname_clean in ['DA','DT','DC','DG','A','T','C','G'] else
                       'rna_base' if resname_clean in ['RA','RU','RC','RG','U'] else
                       'ion'
            }
    
    dna_data = list(residue_dict.values())
    return dna_data

def calculate_quantum_metrics_dna(dna_data):
    """Calculate quantum coherence metrics for DNA systems"""
    if not dna_data:
        return {}
    
    masses = [p['mass'] for p in dna_data]
    n_values = [compute_n_mass_dna(mass) for mass in masses]
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
    superposition_ratio = len([p for p in dna_data if octave_counts[octave_position(compute_n_mass_dna(p['mass']))[2]] > 1]) / len(dna_data)
    
    # Group by molecule type
    type_counts = {}
    dna_bases_present = False
    has_repair_cofactors = False
    
    for item in dna_data:
        mol_type = item['type']
        type_counts[mol_type] = type_counts.get(mol_type, 0) + 1
        if mol_type in ['dna_base', 'rna_base']:
            dna_bases_present = True
        if mol_type == 'cofactor':
            has_repair_cofactors = True
    
    # Quantum information capacity assessment
    base_count = type_counts.get('dna_base', 0) + type_counts.get('rna_base', 0) + type_counts.get('modified_base', 0)
    if base_count >= 4:
        info_capacity = "HIGH"
    elif base_count >= 2:
        info_capacity = "MEDIUM" 
    else:
        info_capacity = "LOW"
    
    return {
        'molecule_count': len(dna_data),
        'type_distribution': type_counts,
        'mass_stats': mass_stats,
        'n_stats': n_stats,
        'octave_distribution': octave_counts,
        'superposition_octaves': superposition_octaves,
        'superposition_ratio': superposition_ratio,
        'mass_tuning': mass_stats['cv'],
        'dna_bases_present': dna_bases_present,
        'base_count': base_count,
        'info_capacity': info_capacity,
        'has_repair_cofactors': has_repair_cofactors,
        'quantum_repair_capable': has_repair_cofactors and superposition_ratio > 0.2
    }

# Global data storage
all_dna_data = []

def load_dna_pdb():
    file_path = filedialog.askopenfilename(filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
    if not file_path:
        return
    dna_molecules = extract_dna_data_pdb(file_path)
    if dna_molecules:
        all_dna_data.extend(dna_molecules)
        update_dna_plot()
        update_dna_analysis(dna_molecules, "PDB")
    else:
        messagebox.showwarning("No Molecules", "No DNA/RNA molecules found.")

def load_dna_cif():
    file_path = filedialog.askopenfilename(filetypes=[("CIF files", "*.cif *.mmcif"), ("All files", "*.*")])
    if not file_path:
        return
    dna_molecules = extract_dna_data_cif(file_path)
    if dna_molecules:
        all_dna_data.extend(dna_molecules)
        update_dna_plot()
        update_dna_analysis(dna_molecules, "CIF")
    else:
        messagebox.showwarning("No Molecules", "No DNA/RNA molecules found.")

def update_dna_analysis(dna_data, file_type):
    analysis_text.delete(1.0, tk.END)
    
    try:
        metrics = calculate_quantum_metrics_dna(dna_data)
        
        text = "DNA/RNA QUANTUM INFORMATION ANALYSIS\n"
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
        
        # Information capacity assessment
        text += f"INFORMATION CAPACITY: {metrics['info_capacity']}\n"
        text += f"  • DNA/RNA bases: {metrics['base_count']}\n"
        if metrics['dna_bases_present']:
            text += "  • Genetic information storage capable ✓\n\n"
        else:
            text += "  • No genetic bases detected\n\n"
        
        # Quantum repair prediction
        if metrics['quantum_repair_capable']:
            text += "QUANTUM REPAIR: LIKELY\n"
            text += "  • Repair cofactors present ✓\n"
            text += "  • Quantum coherence detected ✓\n"
            text += "  • Potential UV damage repair capability\n\n"
        elif metrics['has_repair_cofactors']:
            text += "QUANTUM REPAIR: POSSIBLE\n"
            text += "  • Repair cofactors present ✓\n"
            text += "  • Limited quantum coherence\n\n"
        
        # Base pairing quantum effects
        if metrics['base_count'] >= 2:
            text += "BASE PAIRING QUANTUM EFFECTS:\n"
            text += "  • Multiple bases detected\n"
            text += "  • Potential for quantum information transfer\n"
            text += "  • Possible coherent charge transfer\n\n"
        
        # Spatial organization
        if dna_data:
            n_values = [compute_n_mass_dna(p['mass']) for p in dna_data]
            avg_n = np.mean(n_values)
            spacing, resonance = compute_spatial_radius(avg_n)
            text += f"SPATIAL ORGANIZATION:\n"
            text += f"  Helical spacing: {spacing*1e9:.2f} nm\n"
            text += f"  Resonance factor: {resonance:.3f}\n\n"
        
        n_mean = metrics['n_stats']['mean']
        completed_cycles = int(n_mean)
        current_phase = n_mean - completed_cycles
        
        text += "QUANTUM INFORMATION INTERPRETATION:\n"
        text += f"  Evolutionary depth: ~{completed_cycles} recursive cycles\n"
        text += f"  Current cycle: {current_phase:.1%} complete\n"
        text += f"  Active quantum phases: {len(metrics['octave_distribution'])}/9\n"
        if metrics['superposition_octaves']:
            text += f"  Phase {metrics['superposition_octaves'][0]} shows information coherence\n"
        
        analysis_text.insert(1.0, text)
        
    except Exception as e:
        analysis_text.insert(1.0, f"Error in analysis: {str(e)}")

def update_dna_plot():
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
    
    if all_dna_data:
        color_map = {
            'dna_base': 'blue', 'rna_base': 'red', 
            'modified_base': 'purple', 'backbone': 'green',
            'ion': 'orange', 'cofactor': 'yellow',
            'catalytic': 'brown', 'water': 'cyan'
        }
        
        for molecule in all_dna_data:
            n_val = compute_n_mass_dna(molecule['mass'])
            x, y, octave, angle = octave_position(n_val)
            color = color_map.get(molecule['type'], 'black')
            
            ax.plot(x, y, 'o', color=color, markersize=8, alpha=0.8, markeredgecolor='black', markeredgewidth=1)
            ax.text(x + 0.05, y + 0.05, molecule['resname'], fontsize=7)
    
    if all_dna_data:
        n_values = [compute_n_mass_dna(p['mass']) for p in all_dna_data]
        avg_cycles = np.mean(n_values)
        completed_cycles = int(avg_cycles)
        ax.set_title(f"DNA Quantum Analysis - {completed_cycles} Cycles")
    else:
        ax.set_title("DNA/RNA Quantum Information Analysis")
    
    canvas.draw()

def clear_dna_analysis():
    global all_dna_data
    all_dna_data = []
    analysis_text.delete(1.0, tk.END)
    update_dna_plot()

# Create interface
root = tk.Tk()
root.title("DNA/RNA Quantum Information Analyzer")
root.geometry("1400x800")

# Main container
main_frame = ttk.Frame(root, padding="15")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Control panel
control_frame = ttk.Frame(main_frame, padding="10")
control_frame.grid(row=0, column=0, sticky=(tk.N, tk.S), padx=(0, 10))

ttk.Button(control_frame, text="Load DNA PDB", command=load_dna_pdb, width=15).grid(row=0, column=0, pady=5)
ttk.Button(control_frame, text="Load DNA CIF", command=load_dna_cif, width=15).grid(row=1, column=0, pady=5)
ttk.Button(control_frame, text="Clear", command=clear_dna_analysis, width=15).grid(row=2, column=0, pady=5)

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
update_dna_plot()

root.mainloop()
