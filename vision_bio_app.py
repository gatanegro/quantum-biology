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

# Constants for Bridge Formula
LZ = 1.23488369648610768936904233243207400555883244095073 
QDF = 0.8097928609310675
me = 9.10938356e-31  # kg
pi = np.pi
c = 299792458  # m/s

# Vision system domain references
vision_mass_ref = 3.0e-25  # kg (retinal molecule scale)
vision_energy_ref = 2.3e-19  # J (~1.44 eV - vision photon energy)
vision_radius_ref = 1.2e-9  # m (1.2 nm - photoreceptor spacing)

def F_HQS_tanh(n):
    """Harmonic Quantum Stability - tanh resonance"""
    HQS3 = 0.235543306845346130248413363583896203088885727325439
    return 1.0 + 0.01 * np.tanh(n * HQS3)

def compute_spatial_radius(n, a0_prime=None):
    """Spatial formula for vision systems"""
    if a0_prime is None:
        a0_prime = vision_radius_ref
    
    resonance_factor = F_HQS_tanh(n)
    radius = a0_prime * (LZ ** n) * resonance_factor
    return radius, resonance_factor

def compute_n_mass_vision(mass):
    """Bridge Formula with vision system mass reference"""
    return pi * (np.log10(mass / vision_mass_ref) - np.log10(QDF)) / np.log10(LZ)

def octave_position(n):
    """Map n-value to circular octave position"""
    reduced_n = int(round(n)) % 9
    angle = 2 * pi * reduced_n / 9
    return np.cos(angle), np.sin(angle), reduced_n, angle

def get_mass_by_resname_vision(resname):
    """Vision system molecule database"""
    vision_masses = {
        # Retinal and chromophores
        'RET': 568.87, 'RAL': 568.87, 'RTL': 568.87, 'CBR': 568.87,
        'CRD': 568.87, 'CRO': 568.87,
        # Opsin proteins and key residues
        'OPS': 40000.0, 'OPSD': 40000.0, 'OPSB': 40000.0,  # Approximate
        # Key residues for phototransduction
        'LYS': 146.19,  # Lysine for retinal binding
        'GLU': 147.13,  # Glutamic acid counterion
        'CYS': 121.16,  # Cysteine for disulfide bonds
        'TRP': 204.23,  # Tryptophan for stability
        'TYR': 181.19,  # Tyrosine for hydrogen bonding
        # Lipids and membrane components
        'POPE': 750.0, 'POPC': 750.0, 'CHL': 386.65,
        # Ion channels and signaling
        'CA': 40.08, 'MG': 24.31, 'ZN': 65.38,
        # G-protein components
        'GDP': 443.2, 'GTP': 523.2,
    }
    amu_to_kg = 1.66053906660e-27
    
    # Handle large proteins by using key residue masses
    if resname in ['OPS', 'OPSD', 'OPSB']:
        return 40000.0 * amu_to_kg  # Approximate opsin mass
    else:
        return vision_masses.get(resname, 300.0) * amu_to_kg

def extract_vision_data_pdb(pdb_file):
    """Extract vision-related molecules from PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    vision_data = []
    
    vision_molecules = ['RET', 'RAL', 'RTL', 'CBR', 'CRD', 'CRO', 'OPS', 'OPSD', 'OPSB',
                       'LYS', 'GLU', 'CYS', 'TRP', 'TYR', 'POPE', 'POPC', 'CHL',
                       'CA', 'MG', 'ZN', 'GDP', 'GTP']
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in vision_molecules:
                    mass = get_mass_by_resname_vision(resname)
                    vision_data.append({
                        'resname': resname,
                        'mass': mass,
                        'source': pdb_file,
                        'type': 'chromophore' if resname in ['RET','RAL','RTL','CBR','CRD','CRO'] else
                               'opsin' if resname in ['OPS','OPSD','OPSB'] else
                               'key_residue' if resname in ['LYS','GLU','CYS','TRP','TYR'] else
                               'lipid' if resname in ['POPE','POPC','CHL'] else
                               'ion' if resname in ['CA','MG','ZN'] else
                               'signaling'
                    })
    
    return vision_data

def extract_vision_data_cif(cif_file):
    """Extract vision molecules from CIF file"""
    d = MMCIF2Dict(cif_file)
    vision_data = []
    
    resnames = d.get('_atom_site.label_comp_id', [])
    if isinstance(resnames, str):
        resnames = [resnames]
    
    vision_molecules = ['RET', 'RAL', 'RTL', 'CBR', 'CRD', 'CRO', 'LYS', 'GLU', 
                       'CYS', 'TRP', 'TYR', 'POPE', 'POPC', 'CHL', 'CA', 'MG', 'ZN']
    
    residue_dict = {}
    for resname in resnames:
        resname_clean = resname.strip()
        if resname_clean in vision_molecules and resname_clean not in residue_dict:
            mass = get_mass_by_resname_vision(resname_clean)
            residue_dict[resname_clean] = {
                'resname': resname_clean,
                'mass': mass,
                'source': cif_file,
                'type': 'chromophore' if resname_clean in ['RET','RAL','RTL','CBR','CRD','CRO'] else
                       'key_residue' if resname_clean in ['LYS','GLU','CYS','TRP','TYR'] else
                       'lipid' if resname_clean in ['POPE','POPC','CHL'] else
                       'ion'
            }
    
    vision_data = list(residue_dict.values())
    return vision_data

def calculate_quantum_metrics_vision(vision_data):
    """Calculate quantum coherence metrics for vision systems"""
    if not vision_data:
        return {}
    
    masses = [p['mass'] for p in vision_data]
    n_values = [compute_n_mass_vision(mass) for mass in masses]
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
    superposition_ratio = len([p for p in vision_data if octave_counts[octave_position(compute_n_mass_vision(p['mass']))[2]] > 1]) / len(vision_data)
    
    # Group by molecule type
    type_counts = {}
    chromophore_present = False
    for item in vision_data:
        mol_type = item['type']
        type_counts[mol_type] = type_counts.get(mol_type, 0) + 1
        if mol_type == 'chromophore':
            chromophore_present = True
    
    return {
        'molecule_count': len(vision_data),
        'type_distribution': type_counts,
        'mass_stats': mass_stats,
        'n_stats': n_stats,
        'octave_distribution': octave_counts,
        'superposition_octaves': superposition_octaves,
        'superposition_ratio': superposition_ratio,
        'mass_tuning': mass_stats['cv'],
        'chromophore_present': chromophore_present,
        'vision_capable': chromophore_present and type_counts.get('key_residue', 0) >= 2
    }

# Global data storage
all_vision_data = []

def load_vision_pdb():
    file_path = filedialog.askopenfilename(filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
    if not file_path:
        return
    vision_molecules = extract_vision_data_pdb(file_path)
    if vision_molecules:
        all_vision_data.extend(vision_molecules)
        update_vision_plot()
        update_vision_analysis(vision_molecules, "PDB")
    else:
        messagebox.showwarning("No Molecules", "No vision-related molecules found.")

def load_vision_cif():
    file_path = filedialog.askopenfilename(filetypes=[("CIF files", "*.cif *.mmcif"), ("All files", "*.*")])
    if not file_path:
        return
    vision_molecules = extract_vision_data_cif(file_path)
    if vision_molecules:
        all_vision_data.extend(vision_molecules)
        update_vision_plot()
        update_vision_analysis(vision_molecules, "CIF")
    else:
        messagebox.showwarning("No Molecules", "No vision-related molecules found.")

def update_vision_analysis(vision_data, file_type):
    analysis_text.delete(1.0, tk.END)
    
    try:
        metrics = calculate_quantum_metrics_vision(vision_data)
        
        text = "VISION SYSTEM QUANTUM ANALYSIS\n"
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
        
        # Vision capability prediction
        if metrics['vision_capable']:
            text += "VISION CAPABILITY: HIGH\n"
            text += "  • Chromophore present ✓\n"
            text += "  • Key residues present ✓\n"
            text += "  • Likely functional photoreceptor\n\n"
        elif metrics['chromophore_present']:
            text += "VISION CAPABILITY: MEDIUM\n"
            text += "  • Chromophore present ✓\n"
            text += "  • Limited supporting residues\n\n"
        else:
            text += "VISION CAPABILITY: LOW\n"
            text += "  • No chromophore detected\n\n"
        
        # Photon absorption prediction
        if metrics['chromophore_present']:
            text += "PHOTON INTERACTION:\n"
            text += "  • Retinal chromophore detected\n"
            text += "  • Capable of cis-trans isomerization\n"
            text += "  • Quantum efficiency depends on coherence\n\n"
        
        # Spatial organization
        if vision_data:
            n_values = [compute_n_mass_vision(p['mass']) for p in vision_data]
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

def update_vision_plot():
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
    
    if all_vision_data:
        color_map = {
            'chromophore': 'orange', 'opsin': 'purple', 
            'key_residue': 'blue', 'lipid': 'green', 
            'ion': 'red', 'signaling': 'brown'
        }
        
        for molecule in all_vision_data:
            n_val = compute_n_mass_vision(molecule['mass'])
            x, y, octave, angle = octave_position(n_val)
            color = color_map.get(molecule['type'], 'black')
            
            ax.plot(x, y, 'o', color=color, markersize=8, alpha=0.8, markeredgecolor='black', markeredgewidth=1)
            ax.text(x + 0.05, y + 0.05, molecule['resname'], fontsize=7)
    
    if all_vision_data:
        n_values = [compute_n_mass_vision(p['mass']) for p in all_vision_data]
        avg_cycles = np.mean(n_values)
        completed_cycles = int(avg_cycles)
        ax.set_title(f"Vision System Analysis - {completed_cycles} Cycles")
    else:
        ax.set_title("Vision System Quantum Analysis")
    
    canvas.draw()

def clear_vision_analysis():
    global all_vision_data
    all_vision_data = []
    analysis_text.delete(1.0, tk.END)
    update_vision_plot()

# Create interface
root = tk.Tk()
root.title("Vision System Quantum Analyzer")
root.geometry("1400x800")

# Main container
main_frame = ttk.Frame(root, padding="15")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Control panel
control_frame = ttk.Frame(main_frame, padding="10")
control_frame.grid(row=0, column=0, sticky=(tk.N, tk.S), padx=(0, 10))

ttk.Button(control_frame, text="Load Vision PDB", command=load_vision_pdb, width=15).grid(row=0, column=0, pady=5)
ttk.Button(control_frame, text="Load Vision CIF", command=load_vision_cif, width=15).grid(row=1, column=0, pady=5)
ttk.Button(control_frame, text="Clear", command=clear_vision_analysis, width=15).grid(row=2, column=0, pady=5)

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
update_vision_plot()

root.mainloop()
