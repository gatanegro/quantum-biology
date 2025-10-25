import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import time
import json
from scipy import stats
"""
LOGOS THEORY - Neural Quantum Analyzer
Author: Martin Doina 
"""

# Constants for Bridge Formula
LZ = 1.23488369648610768936904233243207400555883244095073 
QDF = 0.8097928609310675
me = 9.10938356e-31  # kg
pi = np.pi
c = 299792458  # m/s

# BRAIN-SPECIFIC REFERENCES
SYSTEM_REFERENCES = {
    'brain': {
        'mass_ref': 2.8e-25,    # kg (neurotransmitter scale)
        'energy_ref': 3.2e-19,  # J (synaptic transmission range)
        'radius_ref': 2.0e-9,   # m (2 nm - synaptic cleft spacing)
        'description': 'Neural Systems',
        'color': 'purple'
    }
}

# BRAIN MOLECULE DATABASE
MOLECULE_DATABASES = {
    'brain': {
        # Neurotransmitters
        'GLC': 180.16, 'GLU': 147.13, 'GABA': 103.12, 'GLY': 75.07,
        'ASP': 133.10, 'SER': 105.09, 'HIS': 155.16, 'TYR': 181.19,
        'TRP': 204.23, 'PHE': 165.19, 'LYS': 146.19, 'ARG': 174.20,
        
        # Neurotransmitter precursors and metabolites
        'DOPA': 197.19, 'TYR': 181.19, 'TRP': 204.23, 'CHOL': 104.17,
        'SER': 105.09, 'GLN': 146.14, 'ASP': 133.10,
        
        # Receptor components
        'RET': 284.44, 'RAL': 284.44,  # Retinal/chromophores
        'FAD': 785.55, 'FMN': 456.34,  # Flavin cofactors
        'HEM': 616.49, 'HEC': 616.49,  # Heme groups
        
        # Ion channel components
        'CA': 40.08, 'MG': 24.31, 'ZN': 65.38, 'K': 39.10, 'NA': 22.99,
        'CL': 35.45, 'CU': 63.55, 'FE': 55.85, 'MN': 54.94,
        
        # Signaling molecules
        'ATP': 507.18, 'GTP': 523.18, 'CAMP': 329.21, 'CGMP': 345.21,
        'GDP': 443.18, 'GMP': 363.22, 'ADP': 427.20, 'AMP': 347.22,
        
        # Structural proteins
        'TUB': 50000.0, 'ACTN': 100000.0, 'MAP': 80000.0, 'TAU': 45000.0,
        
        # Common residues in neural proteins
        'CYS': 121.16, 'MET': 149.21, 'ASN': 132.12, 'PRO': 115.13,
        'THR': 119.12, 'VAL': 117.15, 'ILE': 131.17, 'LEU': 131.17,
        'ALA': 89.09, 
        
        # Lipids and membrane components
        'CHL': 386.65, 'PLM': 750.0, 'POPE': 750.0, 'POPC': 750.0,
        'PS': 800.0, 'PE': 700.0, 'PC': 700.0,
        
        # Water and ions
        'HOH': 18.02, 'DOD': 20.03, 'SO4': 96.06, 'PO4': 94.97
    }
}

# BRAIN MOLECULE TYPE MAPPINGS
MOLECULE_TYPES = {
    'brain': {
        'neurotransmitter': ['GLC', 'GLU', 'GABA', 'GLY', 'ASP', 'SER', 'HIS', 'TYR', 'TRP', 'PHE', 'LYS', 'ARG', 'DOPA', 'CHOL'],
        'receptor': ['RET', 'RAL', 'FAD', 'FMN', 'HEM', 'HEC'],
        'ion_channel': ['CA', 'MG', 'ZN', 'K', 'NA', 'CL', 'CU', 'FE', 'MN'],
        'signaling': ['ATP', 'GTP', 'CAMP', 'CGMP', 'GDP', 'GMP', 'ADP', 'AMP'],
        'structural': ['TUB', 'ACTN', 'MAP', 'TAU'],
        'membrane': ['CHL', 'PLM', 'POPE', 'POPC', 'PS', 'PE', 'PC'],
        'water_ion': ['HOH', 'DOD', 'SO4', 'PO4']
    }
}

# COLOR SCHEME FOR BRAIN
COLOR_SCHEMES = {
    'brain': {
        'neurotransmitter': 'red',
        'receptor': 'blue', 
        'ion_channel': 'green',
        'signaling': 'orange',
        'structural': 'purple',
        'membrane': 'brown',
        'water_ion': 'cyan'
    }
}

def F_HQS_tanh(n):
    """Harmonic Quantum Stability - tanh resonance"""
    HQS3 = 0.235543306845346130248413363583896203088885727325439
    return 1.0 + 0.01 * np.tanh(n * HQS3)

def compute_spatial_radius(n, system_type):
    """Spatial formula for specific system"""
    ref = SYSTEM_REFERENCES[system_type]
    resonance_factor = F_HQS_tanh(n)
    radius = ref['radius_ref'] * (LZ ** n) * resonance_factor
    return radius, resonance_factor

def compute_n_mass(mass, system_type):
    """Bridge Formula with system-specific mass reference"""
    ref = SYSTEM_REFERENCES[system_type]
    return pi * (np.log10(mass / ref['mass_ref']) - np.log10(QDF)) / np.log10(LZ)

def octave_position(n):
    """Map n-value to circular octave position"""
    reduced_n = int(round(n)) % 9
    angle = 2 * pi * reduced_n / 9
    return np.cos(angle), np.sin(angle), reduced_n, angle

def get_mass_by_resname(resname, system_type):
    """Get mass for molecule in specific system"""
    amu_to_kg = 1.66053906660e-27
    mass = MOLECULE_DATABASES[system_type].get(resname, 300.0)
    return mass * amu_to_kg

def get_molecule_type(resname, system_type):
    """Determine molecule type for specific system"""
    for mol_type, molecules in MOLECULE_TYPES[system_type].items():
        if resname in molecules:
            return mol_type
    return 'water_ion'

def extract_data_pdb(pdb_file, system_type):
    """Extract system-specific molecules from PDB file"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', pdb_file)
        system_data = []
        
        valid_molecules = list(MOLECULE_DATABASES[system_type].keys())
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if resname in valid_molecules:
                        mass = get_mass_by_resname(resname, system_type)
                        mol_type = get_molecule_type(resname, system_type)
                        system_data.append({
                            'resname': resname,
                            'mass': mass,
                            'source': pdb_file,
                            'type': mol_type,
                            'system': system_type
                        })
        
        return system_data
        
    except Exception as e:
        messagebox.showerror("Parsing Error", f"Failed to parse PDB file: {str(e)}")
        return []

def extract_data_cif(cif_file, system_type):
    """Extract system-specific molecules from CIF file"""
    try:
        d = MMCIF2Dict(cif_file)
        system_data = []
        
        resnames = d.get('_atom_site.label_comp_id', [])
        if isinstance(resnames, str):
            resnames = [resnames]
        
        valid_molecules = list(MOLECULE_DATABASES[system_type].keys())
        residue_dict = {}
        
        for resname in resnames:
            resname_clean = resname.strip()
            if resname_clean in valid_molecules and resname_clean not in residue_dict:
                mass = get_mass_by_resname(resname_clean, system_type)
                mol_type = get_molecule_type(resname_clean, system_type)
                residue_dict[resname_clean] = {
                    'resname': resname_clean,
                    'mass': mass,
                    'source': cif_file,
                    'type': mol_type,
                    'system': system_type
                }
        
        system_data = list(residue_dict.values())
        return system_data
        
    except Exception as e:
        messagebox.showerror("Parsing Error", f"Failed to parse CIF file: {str(e)}")
        return []

def calculate_quantum_metrics(system_data, system_type):
    """Calculate quantum coherence metrics for specific system"""
    if not system_data:
        return {}
    
    masses = [p['mass'] for p in system_data]
    n_values = [compute_n_mass(mass, system_type) for mass in masses]
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
    superposition_ratio = len([p for p in system_data if octave_counts[octave_position(compute_n_mass(p['mass'], system_type))[2]] > 1]) / len(system_data)
    
    # Group by molecule type
    type_counts = {}
    for item in system_data:
        mol_type = item['type']
        type_counts[mol_type] = type_counts.get(mol_type, 0) + 1
    
    # System-specific assessments
    system_assessment = assess_system_capability(system_data, system_type, type_counts, superposition_ratio)
    
    return {
        'molecule_count': len(system_data),
        'type_distribution': type_counts,
        'mass_stats': mass_stats,
        'n_stats': n_stats,
        'octave_distribution': octave_counts,
        'superposition_octaves': superposition_octaves,
        'superposition_ratio': superposition_ratio,
        'mass_tuning': mass_stats['cv'],
        'system_assessment': system_assessment,
        'system_type': system_type
    }
def calculate_selectivity_ratio(quantum_metrics):
    # Use the 1.1 nm spacing and 90% coherence
    # Calculate predicted K+:Na+ selectivity
    return predicted_ratio  # Should be ~10,000:1

def calculate_quantum_timescale(superposition_ratio, spacing):
    # Calculate how long quantum states persist
    return coherence_time  # In picoseconds
def calculate_quantum_timescale(superposition_ratio, temperature=310.15):
    """Calculate quantum coherence lifetime in picoseconds"""
    k = 1.380649e-23  # Boltzmann constant
    h = 6.62607015e-34  # Planck constant
    kT = k * temperature  # Thermal energy at body temperature
    
    # Quantum coherence time based on superposition ratio and thermal noise
    base_time = h / (2 * np.pi * kT)  # Fundamental quantum timescale
    coherence_time = base_time * superposition_ratio * 1e12  # Convert to picoseconds
    
    return coherence_time

def calculate_selectivity_ratio(quantum_metrics, ion_type='potassium'):
    """Calculate predicted ion selectivity from quantum parameters"""
    spacing = quantum_metrics.get('optimal_spacing', 2.0e-9)
    superposition = quantum_metrics.get('superposition_ratio', 0.5)
    mass_coherence = quantum_metrics.get('mass_tuning', 20)
    
    # Quantum tunneling probability depends on spacing and coherence
    # Potassium vs sodium selectivity emerges from quantum barrier penetration
    base_selectivity = np.exp((2.0e-9 - spacing) / (0.1e-9))  # Spacing effect
    
    # Quantum coherence enhances discrimination
    coherence_factor = 1.0 + 2.0 * superposition
    
    # Mass tuning indicates evolutionary optimization
    optimization_factor = 1.0 + (mass_coherence / 100.0)
    
    predicted_ratio = base_selectivity * coherence_factor * optimization_factor
    
    if ion_type == 'potassium':
        return int(predicted_ratio * 1000)  # K+:Na+ ratio
    else:
        return int(predicted_ratio)  # General selectivity

def calculate_energy_efficiency(quantum_metrics):
    """Calculate quantum energy savings compared to classical system"""
    superposition = quantum_metrics.get('superposition_ratio', 0.5)
    phase_coherence = len(quantum_metrics.get('superposition_octaves', []))
    mass_tuning = quantum_metrics.get('mass_tuning', 20)
    
    # Quantum systems can process information more efficiently
    # Parallel processing across multiple quantum states reduces energy cost
    base_efficiency = 0.1  # 10% base improvement for quantum systems
    
    # Superposition enables parallel computation
    superposition_gain = superposition * 0.3
    
    # Multiple coherent phases increase computational capacity
    phase_gain = (phase_coherence / 9.0) * 0.4
    
    # Mass tuning indicates evolutionary optimization
    tuning_gain = (mass_tuning / 100.0) * 0.2
    
    total_efficiency = base_efficiency + superposition_gain + phase_gain + tuning_gain
    
    return total_efficiency  # Returns 0.0 to 1.0 (0-100% improvement)

def calculate_biological_predictions(quantum_metrics, system_type):
    """Generate specific biological predictions from quantum analysis"""
    predictions = {}
    
    spacing = quantum_metrics.get('optimal_spacing', 2.0e-9)
    superposition = quantum_metrics.get('superposition_ratio', 0.5)
    coherence_time = calculate_quantum_timescale(superposition)
    
    if system_type == 'brain':
        # Neural-specific predictions
        if spacing < 1.5e-9:
            predictions['ion_selectivity'] = 'HIGH - Optimized for quantum tunneling'
            predictions['selectivity_ratio'] = calculate_selectivity_ratio(quantum_metrics)
        else:
            predictions['ion_selectivity'] = 'MODERATE - Classical filtering dominant'
        
        if coherence_time > 1.0:  # picoseconds
            predictions['timing_precision'] = 'QUANTUM-ENHANCED - Femtosecond precision possible'
        else:
            predictions['timing_precision'] = 'CLASSICAL - Millisecond precision'
        
        if superposition > 0.8:
            predictions['signal_integration'] = 'PARALLEL PROCESSING - Multiple inputs simultaneously'
        else:
            predictions['signal_integration'] = 'SEQUENTIAL PROCESSING - Classical computation'
    
    predictions['energy_efficiency'] = f"{calculate_energy_efficiency(quantum_metrics)*100:.1f}% improvement"
    predictions['coherence_time'] = f"{coherence_time:.1f} ps"
    
    return predictions

def assess_system_capability(system_data, system_type, type_counts, superposition_ratio):
    """Brain-specific capability assessment"""
    assessment = {}
    
    # Neural system assessment
    assessment['neurotransmitters_present'] = type_counts.get('neurotransmitter', 0) > 0
    assessment['receptors_present'] = type_counts.get('receptor', 0) > 0
    assessment['ion_channels_present'] = type_counts.get('ion_channel', 0) > 0
    assessment['signaling_molecules'] = type_counts.get('signaling', 0) > 0
    
    # Determine neural capability level
    essential_components = assessment['neurotransmitters_present'] + assessment['receptors_present'] + assessment['ion_channels_present']
    
    if essential_components >= 3 and superposition_ratio > 0.25:
        assessment['neural_capability'] = 'HIGH'
        assessment['interpretation'] = 'Full synaptic transmission capability'
    elif essential_components >= 2 and superposition_ratio > 0.15:
        assessment['neural_capability'] = 'MEDIUM'
        assessment['interpretation'] = 'Partial neural signaling capability'
    elif essential_components >= 1:
        assessment['neural_capability'] = 'LOW'
        assessment['interpretation'] = 'Basic neural components present'
    else:
        assessment['neural_capability'] = 'MINIMAL'
        assessment['interpretation'] = 'Limited neural function'
    
    # Quantum coherence assessment for neural function
    if superposition_ratio > 0.3:
        assessment['quantum_coherence'] = 'HIGH - Potential for quantum neural processing'
    elif superposition_ratio > 0.2:
        assessment['quantum_coherence'] = 'MEDIUM - Classical neural dynamics dominant'
    else:
        assessment['quantum_coherence'] = 'LOW - Classical signaling only'
    
    return assessment

# Global data storage
all_system_data = []
current_system = 'brain'

def load_pdb():
    global current_system
    file_path = filedialog.askopenfilename(filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")])
    if not file_path:
        return
    
    molecules = extract_data_pdb(file_path, current_system)
    if molecules:
        all_system_data.extend(molecules)
        update_plot()
        update_analysis(molecules, "PDB")
        update_status(f"Loaded {len(molecules)} neural molecules from {file_path}")
    else:
        messagebox.showwarning("No Molecules", f"No neural molecules found in {file_path}")

def load_cif():
    global current_system
    file_path = filedialog.askopenfilename(filetypes=[("CIF files", "*.cif *.mmcif"), ("All files", "*.*")])
    if not file_path:
        return
    
    molecules = extract_data_cif(file_path, current_system)
    if molecules:
        all_system_data.extend(molecules)
        update_plot()
        update_analysis(molecules, "CIF")
        update_status(f"Loaded {len(molecules)} neural molecules from {file_path}")
    else:
        messagebox.showwarning("No Molecules", f"No neural molecules found in {file_path}")

def update_analysis(system_data, file_type):
    analysis_text.delete(1.0, tk.END)
    
    try:
        metrics = calculate_quantum_metrics(system_data, current_system)
        system_info = SYSTEM_REFERENCES[current_system]
        
        text = f"NEURAL QUANTUM ANALYSIS\n"
        text += "=" * 50 + "\n\n"
        
        text += f"STRUCTURE: {file_type}\n"
        text += f"Neural molecules analyzed: {metrics['molecule_count']}\n"
        text += f"Molecule types: {metrics['type_distribution']}\n\n"
        
        text += "NEURAL COMPONENT ANALYSIS:\n"
        assessment = metrics['system_assessment']
        if assessment['neurotransmitters_present']:
            text += "  • Neurotransmitters present ✓\n"
        if assessment['receptors_present']:
            text += "  • Receptor components present ✓\n"
        if assessment['ion_channels_present']:
            text += "  • Ion channel components present ✓\n"
        if assessment['signaling_molecules']:
            text += "  • Signaling molecules present ✓\n"
        text += f"\nNeural Capability: {assessment['neural_capability']}\n"
        text += f"Interpretation: {assessment['interpretation']}\n\n"
        
        text += "QUANTUM NEURAL STATE ANALYSIS:\n"
        text += f"  Recursive depth: {metrics['n_stats']['mean']:.1f} ± {metrics['n_stats']['std']:.1f} cycles\n"
        text += f"  Neural phases: {dict(metrics['octave_distribution'])}\n"
        text += f"  Mass coherence: {metrics['mass_tuning']:.2f}%\n\n"
        
        text += "QUANTUM NEURAL COHERENCE:\n"
        text += f"  Superposition ratio: {metrics['superposition_ratio']:.1%}\n"
        text += f"  Quantum assessment: {assessment['quantum_coherence']}\n"
        if metrics['superposition_octaves']:
            text += f"  Coherent neural phases: {metrics['superposition_octaves']}\n\n"
        else:
            text += "  No neural superposition detected\n\n"
        
        # Spatial organization for neural systems
        if system_data:
                    
            spacing = 2.0e-9  # default
        if system_data:
            n_values = [compute_n_mass(p['mass'], current_system) for p in system_data]
            avg_n = np.mean(n_values)
            spacing, resonance = compute_spatial_radius(avg_n, current_system)
            text += f"SYNAPTIC SPACING ANALYSIS:\n"
            text += f"  Optimal spacing: {spacing*1e9:.1f} nm\n"
            text += f"  Resonance factor: {resonance:.3f}\n\n"
        
        n_mean = metrics['n_stats']['mean']
        completed_cycles = int(n_mean)
        current_phase = n_mean - completed_cycles
        
        text += "QUANTUM NEURAL INTERPRETATION:\n"
        text += f"  Evolutionary depth: ~{completed_cycles} recursive cycles\n"
        text += f"  Current cycle: {current_phase:.1%} complete\n"
        text += f"  Active neural phases: {len(metrics['octave_distribution'])}/9\n"
        if metrics['superposition_octaves']:
            text += f"  Phase {metrics['superposition_octaves'][0]} shows neural quantum coherence\n"
        
        # Additional neural insights
        if assessment['neural_capability'] == 'HIGH':
            text += "\n🎯 This structure shows HIGH potential for quantum-enhanced neural processing!\n"
        elif assessment['neural_capability'] == 'MEDIUM':
            text += "\n⚡ This structure shows classical neural signaling with quantum potential.\n"
        
        analysis_text.insert(1.0, text)
        
    except Exception as e:
        analysis_text.insert(1.0, f"Error in neural analysis: {str(e)}")
        # BIOLOGICAL PREDICTIONS
        text += "\nBIOLOGICAL PREDICTIONS:\n"
        predictions = calculate_biological_predictions({
            'optimal_spacing': spacing,
            'superposition_ratio': metrics['superposition_ratio'],
            'mass_tuning': metrics['mass_tuning'],
            'superposition_octaves': metrics['superposition_octaves']
        }, current_system)
        
        for key, value in predictions.items():
            text += f"  {key.replace('_', ' ').title()}: {value}\n"


def update_plot():
    ax.cla()
    
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.set_aspect('equal')
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    
    # Plot concentric circles
    for radius in [0.25, 0.5, 0.75, 1.0]:
        circle = plt.Circle((0, 0), radius, color='lightgray', fill=False, linestyle='--', alpha=0.5)
        ax.add_artist(circle)
    
    main_circle = plt.Circle((0, 0), 1.0, color='black', fill=False, linewidth=1.5)
    ax.add_artist(main_circle)
    
    # Plot octave lines
    for i in range(9):
        angle = 2 * pi * i / 9
        x = np.cos(angle)
        y = np.sin(angle)
        ax.plot([0, x*1.2], [0, y*1.2], 'gray', linestyle='--', alpha=0.7, linewidth=0.8)
        ax.text(x*1.25, y*1.25, f'O{i}', fontsize=9, color='darkblue')
    
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    
    # Plot neural molecules
    if all_system_data:
        current_data = [m for m in all_system_data if m['system'] == current_system]
        color_scheme = COLOR_SCHEMES[current_system]
        
        for molecule in current_data:
            n_val = compute_n_mass(molecule['mass'], current_system)
            x, y, octave, angle = octave_position(n_val)
            color = color_scheme.get(molecule['type'], 'black')
            
            ax.plot(x, y, 'o', color=color, markersize=8, alpha=0.8, 
                   markeredgecolor='black', markeredgewidth=1)
            ax.text(x + 0.05, y + 0.05, molecule['resname'], fontsize=7)
    
    # Update title
    system_info = SYSTEM_REFERENCES[current_system]
    if all_system_data:
        current_data = [m for m in all_system_data if m['system'] == current_system]
        if current_data:
            n_values = [compute_n_mass(p['mass'], current_system) for p in current_data]
            avg_cycles = np.mean(n_values)
            completed_cycles = int(avg_cycles)
            ax.set_title(f"Neural System - {completed_cycles} Quantum Cycles", 
                        color=system_info['color'], fontsize=14, fontweight='bold')
        else:
            ax.set_title(f"Neural System Analysis - No Data", 
                        color=system_info['color'], fontsize=14, fontweight='bold')
    else:
        ax.set_title(f"Neural Quantum Analysis", 
                    color=system_info['color'], fontsize=14, fontweight='bold')
    
    canvas.draw()

def clear_analysis():
    global all_system_data
    all_system_data = []
    analysis_text.delete(1.0, tk.END)
    update_plot()
    update_status("Neural analysis cleared")

def update_status(message):
    status_var.set(message)

def export_analysis():
    """Export neural analysis results"""
    if not all_system_data:
        messagebox.showwarning("No Data", "No neural data to export.")
        return
    
    try:
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("JSON files", "*.json"), ("All files", "*.*")],
            title="Export Neural Quantum Analysis"
        )
        
        if not file_path:
            return
        
        current_data = [m for m in all_system_data if m['system'] == current_system]
        metrics = calculate_quantum_metrics(current_data, current_system)
        
        if file_path.endswith('.csv'):
            with open(file_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["NEURAL QUANTUM ANALYSIS"])
                writer.writerow(['Export Time', time.strftime('%Y-%m-%d %H:%M:%S')])
                writer.writerow(['System Type', 'Neural/Brain'])
                writer.writerow(['Total Neural Molecules', len(current_data)])
                writer.writerow([])
                
                writer.writerow(['NEURAL ASSESSMENT'])
                assessment = metrics['system_assessment']
                writer.writerow(['Neural Capability', assessment['neural_capability']])
                writer.writerow(['Interpretation', assessment['interpretation']])
                writer.writerow(['Quantum Coherence', assessment['quantum_coherence']])
                writer.writerow([])
                
                writer.writerow(['QUANTUM NEURAL METRICS'])
                writer.writerow(['Superposition Ratio', f"{metrics['superposition_ratio']:.4f}"])
                writer.writerow(['Mass Coherence (%)', f"{metrics['mass_tuning']:.2f}"])
                writer.writerow(['Octave Distribution', str(metrics['octave_distribution'])])
                writer.writerow([])
                
                writer.writerow(['DETAILED NEURAL MOLECULE DATA'])
                writer.writerow(['Resname', 'Type', 'Mass (kg)', 'n-value', 'Octave', 'Neural Role'])
                for molecule in current_data:
                    n_val = compute_n_mass(molecule['mass'], current_system)
                    octave = octave_position(n_val)[2]
                    writer.writerow([
                        molecule['resname'],
                        molecule['type'],
                        f"{molecule['mass']:.2e}",
                        f"{n_val:.4f}",
                        octave,
                        get_neural_role(molecule['resname'])
                    ])
        
        elif file_path.endswith('.json'):
            export_data = {
                'metadata': {
                    'export_time': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'system_type': 'neural_brain',
                    'system_description': 'Neural System Analysis',
                    'neural_molecule_count': len(current_data),
                    'sources': list(set(m['source'] for m in current_data))
                },
                'neural_assessment': metrics['system_assessment'],
                'quantum_metrics': metrics,
                'neural_molecule_data': current_data
            }
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(export_data, f, indent=2, default=str)
        
        messagebox.showinfo("Export Successful", f"Neural data exported to {file_path}")
        update_status(f"Neural analysis exported to {file_path}")
        
    except Exception as e:
        messagebox.showerror("Export Error", f"Failed to export neural data: {str(e)}")

def get_neural_role(resname):
    """Get the neural role for a molecule"""
    roles = {
        'GLC': 'Energy substrate', 'GLU': 'Excitatory neurotransmitter',
        'GABA': 'Inhibitory neurotransmitter', 'GLY': 'Inhibitory neurotransmitter',
        'ASP': 'Excitatory neurotransmitter', 'ATP': 'Energy signaling',
        'GTP': 'G-protein signaling', 'CA': 'Calcium signaling',
        'MG': 'Magnesium block', 'K': 'Potassium channels',
        'NA': 'Sodium channels', 'CL': 'Chloride channels',
        'ZN': 'Synaptic modulation', 'FE': 'Oxygen transport',
        'FAD': 'Redox reactions', 'FMN': 'Redox reactions',
        'RET': 'Visual transduction', 'CHL': 'Membrane integrity'
    }
    return roles.get(resname, 'Neural component')

def export_plot():
    """Export the current neural plot as image"""
    try:
        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), ("All files", "*.*")],
            title="Export Neural Quantum Plot"
        )
        if file_path:
            fig.savefig(file_path, dpi=300, bbox_inches='tight', facecolor='white')
            messagebox.showinfo("Export Successful", f"Neural plot saved to {file_path}")
            update_status(f"Neural plot exported to {file_path}")
    except Exception as e:
        messagebox.showerror("Export Error", f"Failed to export neural plot: {str(e)}")

# Create main interface
root = tk.Tk()
root.title("Neural Quantum Analyzer - LOGOS THEORY")
root.geometry("1400x900")

# Main container
main_frame = ttk.Frame(root, padding="15")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# System selection frame
system_frame = ttk.LabelFrame(main_frame, text="Neural System Analysis", padding="10")
system_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))

# Single system button for brain
ttk.Button(system_frame, text="Neural/Brain System", 
          command=lambda: update_status("Neural system active"),
          width=25, style='Accent.TButton').grid(row=0, column=0, padx=5)

# Control panel
control_frame = ttk.Frame(main_frame, padding="10")
control_frame.grid(row=1, column=0, sticky=(tk.N, tk.S), padx=(0, 10))

# Control buttons
ttk.Button(control_frame, text="Load Neural PDB", command=load_pdb, width=18).grid(row=0, column=0, pady=5)
ttk.Button(control_frame, text="Load Neural CIF", command=load_cif, width=18).grid(row=1, column=0, pady=5)
ttk.Button(control_frame, text="Clear Analysis", command=clear_analysis, width=18).grid(row=2, column=0, pady=5)
ttk.Button(control_frame, text="Export Neural Data", command=export_analysis, width=18).grid(row=3, column=0, pady=5)
ttk.Button(control_frame, text="Export Neural Plot", command=export_plot, width=18).grid(row=4, column=0, pady=5)

# Analysis display
analysis_frame = ttk.LabelFrame(control_frame, text="Neural Quantum Analysis", padding="5")
analysis_frame.grid(row=5, column=0, pady=10, sticky=(tk.W, tk.E))

analysis_text = tk.Text(analysis_frame, width=35, height=25, wrap=tk.WORD, 
                       font=('Courier', 9), bg='white')
analysis_text.pack(fill=tk.BOTH, expand=True)

# Plot area
plot_frame = ttk.Frame(main_frame)
plot_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))

# Create figure
fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
fig.patch.set_facecolor('white')

canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# Status bar
status_var = tk.StringVar(value="Ready - Load a neural structure file (PDB/CIF) to begin analysis")
status_bar = ttk.Label(root, textvariable=status_var, relief=tk.SUNKEN, anchor=tk.W)
status_bar.grid(row=1, column=0, sticky=(tk.W, tk.E))

# Configure grid weights
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
main_frame.columnconfigure(1, weight=1)
main_frame.rowconfigure(1, weight=1)
plot_frame.columnconfigure(0, weight=1)
plot_frame.rowconfigure(0, weight=1)

# Initialize
analysis_text.insert(1.0, "NEURAL QUANTUM ANALYZER\n\n" +
                    "• Specifically designed for brain protein analysis\n" +
                    "• Analyzes neurotransmitters, receptors, ion channels\n" +
                    "• Assesses quantum coherence in neural systems\n\n" +
                    "SUGGESTED PDB FILES:\n" +
                    "• NMDA/AMPA receptors\n" +
                    "• GABA receptors\n" +
                    "• Ion channels\n" +
                    "• Neurotransmitter transporters\n\n" +
                    "Load a neural structure to begin!")

update_plot()

root.mainloop()
