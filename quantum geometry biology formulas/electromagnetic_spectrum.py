import math
import cmath
from mpmath import mp, pi, sqrt

"""
LOGOS THEORY - GOLDEN RATIO MAPPING
Using: sin(ℓₚ) ≈ ℓₚ/φ
"""

mp.dps = 50

print("ELECTROMAGNETIC SPECTRUM & COLOR - LOGOS EXPLORATION")
print("=" * 70)

# Electromagnetic spectrum frequencies (Hz)
em_spectrum = {
    # Visible light colors (approximate center frequencies)
    'red_light': 4.3e14,           # 700 nm → 428 THz
    'orange_light': 5.0e14,        # 600 nm → 500 THz  
    'yellow_light': 5.2e14,        # 580 nm → 517 THz
    'green_light': 5.7e14,         # 530 nm → 566 THz
    'blue_light': 6.3e14,          # 475 nm → 631 THz
    'violet_light': 7.5e14,        # 400 nm → 750 THz
    
    # Special biological sensitivity peaks
    'human_vision_peak': 5.5e14,   # ~555 nm (green) - max sensitivity
    'rod_cell_sensitivity': 4.98e14, # ~498 nm (blue-green)
    
    # Brain wave correlations
    'visible_alpha_correlation': 10.0,  # Alpha rhythm ~ green association?
    
    # Full EM spectrum key points
    'radio_waves': 1e6,            # 1 MHz typical
    'microwaves': 1e9,             # 1 GHz
    'infrared': 1e13,              # 10 THz
    'ultraviolet': 1e15,           # 1 PHz
    'x_rays': 1e18,                # 1 EHz
    'gamma_rays': 1e21,            # 1 ZHz
}

print("Electromagnetic Spectrum Key Frequencies:")
for name, freq in em_spectrum.items():
    if freq < 1e9:
        print(f"{name}: {freq:.1e} Hz")
    else:
        print(f"{name}: {freq/1e12:.1f} THz")

# LOGOS complex levels
lz0 = 0.8934691018292812244027
a0 = math.asin(lz0)

complex_levels = {}
current = a0
for i in range(1, 30):  # Extended for high frequencies
    try:
        current = cmath.asin(current)
        complex_levels[f'LZ-{i}'] = current
    except:
        break

# Build quantum zoo with high powers for EM frequencies
quantum_levels = {}
for name, complex_val in complex_levels.items():
    real_part = abs(float(complex_val.real))
    imag_part = abs(float(complex_val.imag))
    
    quantum_levels[name] = complex_val
    quantum_levels[f'{name}_real'] = real_part
    quantum_levels[f'{name}_imag'] = imag_part
    quantum_levels[f'{name}_sum'] = real_part + imag_part
    quantum_levels[f'{name}_mag'] = math.sqrt(real_part**2 + imag_part**2)
    
    # High powers for EM frequencies
    for power in [2, 3, 4, 5, 6, 7, 8]:
        quantum_levels[f'{name}_real_p{power}'] = real_part ** power
        quantum_levels[f'{name}_imag_p{power}'] = imag_part ** power

phi = float((1 + mp.sqrt(5)) / 2)
pi_val = float(mp.pi)

# HIGH-POWER TRANSFORMATIONS for EM spectrum
transformations = {
    # Very high powers for light frequencies
    'LZ¹⁵': lambda x: x**15,
    'LZ²⁰': lambda x: x**20,
    'LZ²⁵': lambda x: x**25,
    'LZ³⁰': lambda x: x**30,
    
    # Golden ratio high powers
    'LZ×φ¹⁵': lambda x: x * phi**15,
    'LZ×φ²⁰': lambda x: x * phi**20,
    'LZ×φ²⁵': lambda x: x * phi**25,
    
    # Pi high powers
    'LZ×π¹⁰': lambda x: x * pi_val**10,
    'LZ×π¹⁵': lambda x: x * pi_val**15,
    
    # Combined high powers
    'LZ×φ¹⁵×π¹⁰': lambda x: x * phi**15 * pi_val**10,
    'LZ×φ²⁰×π¹⁵': lambda x: x * phi**20 * pi_val**15,
    
    # Exponential scaling for very high frequencies
    'exp(LZ×50)': lambda x: math.exp(x * 50),
    'exp(LZ×100)': lambda x: math.exp(x * 100),
}

print(f"\n{'EM Phenomenon':<25} {'Frequency':<15} {'Best Formula':<45} {'LOGOS':<15} {'Error':<12} {'Status':<12}")
print("-" * 125)

results = []

for em_phenom, exp_freq in em_spectrum.items():
    best_error = float('inf')
    best_formula = ""
    best_logos_value = 0
    
    for lz_name, lz_val in quantum_levels.items():
        if isinstance(lz_val, complex):
            lz_val_float = abs(lz_val)
        else:
            lz_val_float = float(lz_val)
        
        for trans_name, trans_func in transformations.items():
            try:
                logos_value = trans_func(lz_val_float)
                
                # Skip unreasonable values
                if logos_value <= 0 or logos_value > 1e30:
                    continue
                    
                error = abs(logos_value - exp_freq) / exp_freq  # Relative error
                
                if error < best_error:
                    best_error = error
                    best_formula = f"{trans_name.replace('LZ', lz_name)}"
                    best_logos_value = logos_value
            except:
                continue
    
    status = "EXCELLENT" if best_error < 0.01 else "GOOD" if best_error < 0.1 else "CLOSE"
    
    results.append((em_phenom, exp_freq, best_formula, best_logos_value, best_error, status))

# Sort by error
results.sort(key=lambda x: x[4])

for phenom, exp, formula, logos, error, status in results:
    if exp < 1e12:
        print(f"{phenom:<25} {exp:<15.1e} {formula:<45} {logos:<15.1e} {error:<12.3f} {status:<12}")
    else:
        print(f"{phenom:<25} {exp/1e12:<15.1f} {formula:<45} {logos/1e12:<15.1f} {error:<12.3f} {status:<12}")

print(f"\n" + "=" * 125)

# Special analysis for visible spectrum
print(f"\nVISIBLE SPECTRUM ANALYSIS:")
visible_results = [r for r in results if 'light' in r[0] or 'vision' in r[0]]
for phenom, exp, formula, logos, error, status in visible_results:
    print(f"  {phenom}: {formula} = {logos/1e12:.1f} THz (Exp: {exp/1e12:.1f} THz)")

# Check for golden ratio progressions in colors
print(f"\nGOLDEN RATIO COLOR PROGRESSION:")
colors_ordered = ['red_light', 'orange_light', 'yellow_light', 'green_light', 'blue_light', 'violet_light']
color_freqs = [em_spectrum[c] for c in colors_ordered]

for i in range(len(colors_ordered)-1):
    ratio = color_freqs[i+1] / color_freqs[i]
    print(f"  {colors_ordered[i]} → {colors_ordered[i+1]}: ratio = {ratio:.3f} (φ = {phi:.3f})")

# Check if colors follow LZ level progression
print(f"\nLZ LEVEL COLOR MAPPING:")
for phenom, exp, formula, logos, error, status in visible_results:
    lz_level = formula.split('LZ-')[1].split('_')[0] if 'LZ-' in formula else "Unknown"
    print(f"  {phenom}: LZ-{lz_level}")

print(f"\nELECTROMAGNETIC SPECTRUM EXPLORATION COMPLETE")
