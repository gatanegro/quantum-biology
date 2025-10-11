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

# Biological domain references
bio_mass_ref = 8.9e-25  # kg (chlorophyll-scale)
bio_energy_ref = 2.88e-19  # J (1.8 eV)
bio_radius_ref = 1.5e-9  # m (1.5 nm)

def F_HQS_tanh(n):
    """Harmonic Quantum Stability - tanh resonance"""
    HQS3 = 0.235543306845346130248413363583896203088885727325439
    return 1.0 + 0.01 * np.tanh(n * HQS3)

def F_HQS_sin(n):
    """Harmonic Quantum Stability - sine resonance"""
    HQS3 = 0.235543306845346130248413363583896203088885727325439
    return 1.0 + 0.005 * np.sin(n * HQS3)

def F_HQS_exp(n):
    """Harmonic Quantum Stability - exponential resonance"""
    HQS3 = 0.235543306845346130248413363583896203088885727325439
    return 1.0 + 0.02 * np.exp(-abs(n * HQS3))

def compute_spatial_radius_all(n, a0_prime=None, resonance_type='tanh'):
    """Complete spatial formula with resonance options"""
    if a0_prime is None:
        a0_prime = bio_radius_ref
    
    if resonance_type == 'tanh':
        resonance_factor = F_HQS_tanh(n)
    elif resonance_type == 'sin':
        resonance_factor = F_HQS_sin(n)
    elif resonance_type == 'exp':
        resonance_factor = F_HQS_exp(n)
    else:
        resonance_factor = 1.0
    
    radius = a0_prime * (LZ ** n) * resonance_factor
    return radius, resonance_factor

def compute_n_mass_bio(mass):
    """Bridge Formula with biological mass reference"""
    return pi * (np.log10(mass / bio_mass_ref) - np.log10(QDF)) / np.log10(LZ)

def compute_n_energy_bio(energy):
    """Bridge Formula with biological energy reference""" 
    return pi * (np.log10(energy / bio_energy_ref) - np.log10(QDF)) / np.log10(LZ)

def octave_position(n):
    """Map n-value to circular octave position"""
    reduced_n = int(round(n)) % 9
    angle = 2 * pi * reduced_n / 9
    return np.cos(angle), np.sin(angle), reduced_n, angle

def get_mass_by_resname(resname):
    """Pigment mass database"""
    pigment_masses = {
        'CLA': 893.49, 'CHL': 907.47, 'CHLA': 893.49, 
        'BCL': 907.47, 'BCLA': 907.47, 'CAR': 536.87, 
        'BCR': 560.0, 'BCT': 536.87, 'LUT': 568.87, 
        'ZEA': 568.87, 'VAU': 600.0, 'NEX': 600.0,
        'PL9': 750.0, 'PQN': 750.0, 'PHO': 800.0,
        'UQ1': 750.0, 'UDQ': 750.0, 'HEM': 616.49, 
        'HEC': 616.49, 'LMG': 750.0, 'DGK': 800.0, 
        'CRD': 750.0,
    }
    amu_to_kg = 1.66053906660e-27
    return pigment_masses.get(resname, 500.0) * amu_to_kg

def extract_pigment_data_pdb(pdb_file):
    """Extract pigments from PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    pigment_data = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                mass = get_mass_by_resname(resname)
                if mass != 500.0 * 1.66053906660e-27:
                    pigment_data.append({
                        'resname': resname,
                        'mass': mass,
                        'source': pdb_file
                    })
    return pigment_data

def extract_pigment_data_cif(cif_file):
    """Extract pigments from CIF file"""
    d = MMCIF2Dict(cif_file)
    pigment_data = []
    
    resnames = d.get('_atom_site.label_comp_id', [])
    if isinstance(resnames, str):
        resnames = [resnames]
    
    residue_dict = {}
    for resname in resnames:
        resname_clean = resname.strip()
        if resname_clean not in residue_dict:
            mass = get_mass_by_resname(resname_clean)
            if mass != 500.0 * 1.66053906660e-27:
                residue_dict[resname_clean] = {
                    'resname': resname_clean,
                    'mass': mass,
                    'source': cif_file
                }
    
    pigment_data = list(residue_dict.values())
    return pigment_data

def calculate_quantum_metrics(pigment_data):
    """Calculate quantum coherence metrics"""
    if not pigment_data:
        return {}
    
    masses = [p['mass'] for p in pigment_data]
    n_values = [compute_n_mass_bio(mass) for mass in masses]
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
    superposition_ratio = len([p for p in pigment_data if octave_counts[octave_position(compute_n_mass_bio(p['mass']))[2]] > 1]) / len(pigment_data)
    
    if len(set(octaves)) > 1:
        entropy = stats.entropy([octave_counts.get(i, 0) for i in range(9)])
    else:
        entropy = 0
    
    return {
        'pigment_count': len(pigment_data),
        'mass_stats': mass_stats,
        'n_stats': n_stats,
        'octave_distribution': octave_counts,
        'superposition_octaves': superposition_octaves,
        'superposition_ratio': superposition_ratio,
        'entropy': entropy,
        'quantum_coherence': superposition_ratio * 100,
        'mass_tuning': mass_stats['cv']
    }

def analyze_complete_spatial_unification(pigment_data):
    """Complete analysis with all resonance function options"""
    results = []
    
    experimental_energies = {
        'CLA': 1.77, 'CHL': 1.77, 'CHLA': 1.77,
        'CAR': 2.38, 'BCR': 2.38, 'BCT': 2.38,
        'PHO': 1.78, 'LUT': 2.42, 'ZEA': 2.42, 
        'VAU': 2.35, 'HEM': 2.30, 'PL9': 2.10,
    }
    
    for pigment in pigment_data:
        mass = pigment['mass']
        resname = pigment['resname']
        
        n_mass = compute_n_mass_bio(mass)
        octave_mass = octave_position(n_mass)[2]
        
        radius_tanh, resonance_tanh = compute_spatial_radius_all(n_mass, resonance_type='tanh')
        radius_sin, resonance_sin = compute_spatial_radius_all(n_mass, resonance_type='sin')
        radius_exp, resonance_exp = compute_spatial_radius_all(n_mass, resonance_type='exp')
        
        energy_analysis = {}
        if resname in experimental_energies:
            energy_eV = experimental_energies[resname]
            energy_J = energy_eV * 1.602e-19
            n_energy = compute_n_energy_bio(energy_J)
            octave_energy = octave_position(n_energy)[2]
            energy_analysis = {
                'energy_eV': energy_eV,
                'n_energy': n_energy,
                'octave_energy': octave_energy,
                'mass_energy_alignment': abs(n_mass - n_energy) < 2.0
            }
        
        results.append({
            'resname': resname,
            'mass_kg': mass,
            'n_mass': n_mass,
            'octave_mass': octave_mass,
            'optimal_radius_tanh_m': radius_tanh,
            'optimal_radius_sin_m': radius_sin,
            'optimal_radius_exp_m': radius_exp,
            'resonance_tanh': resonance_tanh,
            'resonance_sin': resonance_sin,
            'resonance_exp': resonance_exp,
            'energy_analysis': energy_analysis,
        })
    
    return results

def calculate_spatial_optimization_metrics_all(results):
    """Calculate spatial optimization metrics for all resonance types"""
    if not results or len(results) < 2:
        return {
            'tanh': {'avg_radius_nm': 0, 'spatial_optimization': 0, 'resonance_stability': 0, 'radius_std_nm': 0},
            'sin': {'avg_radius_nm': 0, 'spatial_optimization': 0, 'resonance_stability': 0, 'radius_std_nm': 0},
            'exp': {'avg_radius_nm': 0, 'spatial_optimization': 0, 'resonance_stability': 0, 'radius_std_nm': 0}
        }
    
    radii_tanh = [r['optimal_radius_tanh_m'] for r in results]
    radii_sin = [r['optimal_radius_sin_m'] for r in results]
    radii_exp = [r['optimal_radius_exp_m'] for r in results]
    
    reasonable_tanh = sum(1 for r in radii_tanh if 0.5e-9 <= r <= 5e-9)
    reasonable_sin = sum(1 for r in radii_sin if 0.5e-9 <= r <= 5e-9)
    reasonable_exp = sum(1 for r in radii_exp if 0.5e-9 <= r <= 5e-9)
    
    # Handle single pigment case
    if len(results) == 1:
        stability_tanh = 1.0
        stability_sin = 1.0
        stability_exp = 1.0
        radius_std_tanh = 0
        radius_std_sin = 0
        radius_std_exp = 0
    else:
        resonance_tanh = [r['resonance_tanh'] for r in results]
        resonance_sin = [r['resonance_sin'] for r in results]
        resonance_exp = [r['resonance_exp'] for r in results]
        
        stability_tanh = 1.0 - np.std(resonance_tanh) if len(resonance_tanh) > 1 else 1.0
        stability_sin = 1.0 - np.std(resonance_sin) if len(resonance_sin) > 1 else 1.0
        stability_exp = 1.0 - np.std(resonance_exp) if len(resonance_exp) > 1 else 1.0
        
        radius_std_tanh = np.std(radii_tanh) if len(radii_tanh) > 1 else 0
        radius_std_sin = np.std(radii_sin) if len(radii_sin) > 1 else 0
        radius_std_exp = np.std(radii_exp) if len(radii_exp) > 1 else 0
    
    return {
        'tanh': {
            'avg_radius_nm': np.mean(radii_tanh) * 1e9,
            'spatial_optimization': reasonable_tanh / len(results),
            'resonance_stability': max(0, stability_tanh),
            'radius_std_nm': radius_std_tanh * 1e9
        },
        'sin': {
            'avg_radius_nm': np.mean(radii_sin) * 1e9,
            'spatial_optimization': reasonable_sin / len(results),
            'resonance_stability': max(0, stability_sin),
            'radius_std_nm': radius_std_sin * 1e9
        },
        'exp': {
            'avg_radius_nm': np.mean(radii_exp) * 1e9,
            'spatial_optimization': reasonable_exp / len(results),
            'resonance_stability': max(0, stability_exp),
            'radius_std_nm': radius_std_exp * 1e9
        }
    }

# Global data storage
all_pigment_data = []

def load_pdb_and_plot():
    file_path = filedialog.askopenfilename(filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
    if not file_path:
        return
    pigments = extract_pigment_data_pdb(file_path)
    if pigments:
        all_pigment_data.extend(pigments)
        update_plot()
        update_analysis_display(pigments, "PDB")
    else:
        messagebox.showwarning("No Pigments", "No photosynthetic pigments found.")

def load_cif_and_plot():
    file_path = filedialog.askopenfilename(filetypes=[("CIF files", "*.cif *.mmcif"), ("All files", "*.*")])
    if not file_path:
        return
    pigments = extract_pigment_data_cif(file_path)
    if pigments:
        all_pigment_data.extend(pigments)
        update_plot()
        update_analysis_display(pigments, "CIF")
    else:
        messagebox.showwarning("No Pigments", "No photosynthetic pigments found.")

def update_analysis_display(pigments, file_type):
    analysis_text.delete(1.0, tk.END)
    
    try:
        metrics = calculate_quantum_metrics(pigments)
        
        text = "SPIRAL RECURSIVE UNIFICATION\n"
        text += "=" * 50 + "\n\n"
        
        text += f"STRUCTURE: {file_type}\n"
        text += f"Pigments analyzed: {metrics['pigment_count']}\n\n"
        
        text += "SPIRAL RECURSION STATE:\n"
        text += f"  Recursive depth: {metrics['n_stats']['mean']:.1f} ± {metrics['n_stats']['std']:.1f} cycles\n"
        text += f"  Current cycle phases: {dict(metrics['octave_distribution'])}\n"
        text += f"  Mass coherence: {metrics['mass_tuning']:.2f}%\n\n"
        
        text += "QUANTUM COHERENCE:\n"
        text += f"  Superposition ratio: {metrics['superposition_ratio']:.1%}\n"
        if metrics['superposition_octaves']:
            text += f"  Coherent phases: {metrics['superposition_octaves']}\n\n"
        else:
            text += "  No superposition detected\n\n"
        
        # Simple spatial calculation
        if pigments:
            n_values = [compute_n_mass_bio(p['mass']) for p in pigments]
            avg_n = np.mean(n_values)
            simple_spacing = compute_spatial_radius_all(avg_n)[0] * 1e9
            text += f"SPATIAL ORGANIZATION:\n"
            text += f"  Average spacing: {simple_spacing:.1f} nm\n\n"
        
        n_mean = metrics['n_stats']['mean']
        completed_cycles = int(n_mean)
        current_phase = n_mean - completed_cycles
        
        text += "INTERPRETATION:\n"
        text += f"  Evolutionary depth: ~{completed_cycles} recursive cycles\n"
        text += f"  Current cycle: {current_phase:.1%} complete\n"
        text += f"  Active phases: {len(metrics['octave_distribution'])}/9\n"
        if metrics['superposition_octaves']:  # FIX: Check if list is not empty
            text += f"  Phase {metrics['superposition_octaves'][0]} dominant coherence\n"
        else:
            text += "  No dominant coherence phase\n"
        
        analysis_text.insert(1.0, text)
        
    except Exception as e:
        analysis_text.insert(1.0, f"Error: {str(e)}\n\nTry loading a different structure.")

def update_plot():
    ax.cla()
    
    # Simple styling - no font issues
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
    
    if all_pigment_data:
        color_map = {
            'CLA': 'green', 'CHL': 'darkgreen', 'CHLA': 'lime',
            'CAR': 'orange', 'BCR': 'red', 'BCT': 'darkred',
            'PL9': 'blue', 'PHO': 'purple', 'HEM': 'brown',
            'LUT': 'yellow', 'ZEA': 'gold'
        }
        
        for pigment in all_pigment_data:
            n_val = compute_n_mass_bio(pigment['mass'])
            x, y, octave, angle = octave_position(n_val)
            color = color_map.get(pigment['resname'], 'black')
            
            ax.plot(x, y, 'o', color=color, markersize=8, alpha=0.8, markeredgecolor='black', markeredgewidth=1)
            ax.text(x + 0.05, y + 0.05, pigment['resname'], fontsize=7)
    
    if all_pigment_data:
        n_values = [compute_n_mass_bio(p['mass']) for p in all_pigment_data]
        avg_cycles = np.mean(n_values)
        completed_cycles = int(avg_cycles)
        ax.set_title(f"Spiral Analysis - {completed_cycles} Cycles")
    else:
        ax.set_title("Spiral Analysis")
    
    canvas.draw()

def clear_analysis():
    global all_pigment_data
    all_pigment_data = []
    analysis_text.delete(1.0, tk.END)
    update_plot()

# Create interface
root = tk.Tk()
root.title("Spiral Recursive Analysis")
root.geometry("1400x800")

# Main container
main_frame = ttk.Frame(root, padding="15")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Control panel
control_frame = ttk.Frame(main_frame, padding="10")
control_frame.grid(row=0, column=0, sticky=(tk.N, tk.S), padx=(0, 10))

# Simple buttons - no font issues
ttk.Button(control_frame, text="Load PDB", command=load_pdb_and_plot, width=15).grid(row=0, column=0, pady=5)
ttk.Button(control_frame, text="Load CIF", command=load_cif_and_plot, width=15).grid(row=1, column=0, pady=5)
ttk.Button(control_frame, text="Clear", command=clear_analysis, width=15).grid(row=2, column=0, pady=5)

# Analysis display
analysis_text = tk.Text(control_frame, width=35, height=25, wrap=tk.WORD, font=('Courier', 9), bg='white')
analysis_text.grid(row=3, column=0, pady=10)

# Plot area
plot_frame = ttk.Frame(main_frame)
plot_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))

# Create figure with default fonts
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
update_plot()

root.mainloop()
