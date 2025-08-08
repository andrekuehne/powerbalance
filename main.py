# %% Power Balance Example - Python Version
import sys
sys.path.append('src')  # add src directory to path for module imports
from PowerBalance import PowerBalanceCalculator
import h5py
import numpy as np
import skrf as rf
import matplotlib.pyplot as plt

# %% Load data from HDF5 files
print("Loading field data...")
with h5py.File('example/2ch_loops.exported/e_field.mat', 'r') as f:
    print("E-field keys:", list(f.keys()))
    tmp = f['e'][:]
    e = tmp['real'] + 1j*tmp['imag']  # convert to complex
    e = e.transpose(4, 3, 2, 1, 0).astype(np.complex128)  # Reverse all dimensions

with h5py.File('example/2ch_loops.exported/h_field.mat', 'r') as f:
    print("H-field keys:", list(f.keys()))
    tmp = f['h'][:]
    h = tmp['real'] + 1j*tmp['imag']  # convert to complex
    h = h.transpose(4, 3, 2, 1, 0).astype(np.complex128)  # Reverse all dimensions

with h5py.File('example/2ch_loops.exported/current_density.mat', 'r') as f:
    print("Current density keys:", list(f.keys()))
    tmp = f['c'][:]
    c = tmp['real'] + 1j*tmp['imag']  # convert to complex
    c = c.transpose(4, 3, 2, 1, 0).astype(np.complex128)  # Reverse all dimensions

with h5py.File('example/2ch_loops.exported/u_matrix.mat', 'r') as f:
    print("U-matrix keys:", list(f.keys()))
    tmp = f['uMat'][:]
    u_matrix = tmp['real'] + 1j*tmp['imag']  # convert to complex
    u_matrix = u_matrix.T.astype(np.complex128)

with h5py.File('example/2ch_loops.exported/i_matrix.mat', 'r') as f:
    print("I-matrix keys:", list(f.keys()))
    tmp = f['iMat'][:]
    i_matrix = tmp['real'] + 1j*tmp['imag']  # convert to complex
    i_matrix = i_matrix.T.astype(np.complex128)


with h5py.File('example/2ch_loops.exported/meshdata.mat', 'r') as f:
    print("Mesh data keys:", list(f.keys()))
    x_line = f['xLine'][:].ravel().astype(np.float64)  # ensure 1D
    y_line = f['yLine'][:].ravel().astype(np.float64)  # ensure 1D  
    z_line = f['zLine'][:].ravel().astype(np.float64)  # ensure 1D

print(f"Original field shapes: e={e.shape}, h={h.shape}, c={c.shape}")

# %% Crop last points in field data (match MATLAB indexing)
e = e[:, :-1, :-1, :-1, :]
h = h[:, :-1, :-1, :-1, :]
c = c[:, :-1, :-1, :-1, :]

print(f"Cropped field shapes: e={e.shape}, h={h.shape}, c={c.shape}")

# %% Load S-parameters
rf_net = rf.network.Network('example/2ch_loops.exported/s_matrix.s2p')
s_matrix = rf_net.s[0, :, :]  # extract S-parameter matrix at first frequency

print(f"S-matrix shape: {s_matrix.shape}")

# %% Calculate conductivity from current density and electric field
# Use only first port for conductivity calculation (as in MATLAB)
sigma = c[0] / e[0]  # element-wise division
sigma = np.real(sigma)  # take real part
sigma[~np.isfinite(sigma)] = 0  # replace NaN/inf with 0

print(f"Sigma shape: {sigma.shape}")
print(f"Sigma range: {np.min(sigma)} to {np.max(sigma)}")

# %% Generate binary masks for material differentiation
mask_fr4 = (sigma > 0) & (sigma < 0.001)
mask_phantom = (sigma > 0) & (sigma >= 0.001)

sigma_masks = [mask_fr4, mask_phantom]

print(f"FR4 mask has {np.sum(mask_fr4)} active cells")
print(f"Phantom mask has {np.sum(mask_phantom)} active cells")

# %% Generate binary masks for lumped component differentiation
# [1 1 0 0 0 0 0 0 0] - matching caps
# [0 0 1 1 1 1 1 1 0] - fixed & tuning caps  
# [0 0 0 0 0 0 0 0 1] - decoupling cap
mask_cm = np.array([True, True, False, False, False, False, False, False, False])  # matching caps
mask_ct = np.array([False, False, True, True, True, True, True, True, False])     # fixed & tuning caps
mask_cd = np.array([False, False, False, False, False, False, False, False, True]) # decoupling cap

lumped_masks = [mask_cm, mask_ct, mask_cd]

print(f"Lumped element masks created for {u_matrix.shape[0]} elements")

# %% Calculate power balance
print("Calculating power balance...")
forward_power = 0.02  # 20 mW forward power
far_field_offset = 2

calculator = PowerBalanceCalculator(far_field_offset=far_field_offset)

grid_vectors = (x_line.astype(np.float32), 
                y_line.astype(np.float32), 
                z_line.astype(np.float32))

result = calculator.calculate_power_balance(
    e, h, sigma, grid_vectors,
    u_matrix, i_matrix, s_matrix, forward_power,
    sigma_masks=sigma_masks, lumped_masks=lumped_masks
)

print("Power balance calculation complete!")
print(f"Available matrices: fwd, cpl, lmp, mat, rad")
print(f"Material masks: {len(result.mat_masked) if result.mat_masked else 0}")
print(f"Lumped masks: {len(result.lmp_masked) if result.lmp_masked else 0}")

# %% Example plot - Power balance under phase-only shimming
print("Generating power balance plot...")

n_points = 361
phases = np.linspace(0, 180, n_points)  # phase sweep from 0 to 180 degrees

# Create excitation vectors (as in MATLAB)
# v = [1; exp(1j*phase)] / sqrt(2)
v_vectors = np.zeros((2, n_points), dtype=complex)
v_vectors[0, :] = 1  # first port always 1
v_vectors[1, :] = np.exp(1j * phases * np.pi / 180)  # second port with phase
v_vectors = v_vectors / np.sqrt(2)  # normalize

# Initialize power arrays
p_ref = np.zeros(n_points)
p_rad = np.zeros(n_points)
p_phantom = np.zeros(n_points)
p_substrate = np.zeros(n_points)
p_match = np.zeros(n_points)
p_fixed = np.zeros(n_points)
p_decouple = np.zeros(n_points)

# Calculate power for each phase point
print("Computing power vs phase...")
for k in range(n_points):
    v = v_vectors[:, k]
    
    # Calculate quadratic forms: real(v^H * Matrix * v)
    p_ref[k] = np.real(v.conj().T @ result.cpl @ v)
    p_rad[k] = np.real(v.conj().T @ result.rad @ v)
    
    if result.mat_masked:
        p_substrate[k] = np.real(v.conj().T @ result.mat_masked[0] @ v)  # FR4
        p_phantom[k] = np.real(v.conj().T @ result.mat_masked[1] @ v)    # Phantom
    
    if result.lmp_masked:
        p_match[k] = np.real(v.conj().T @ result.lmp_masked[0] @ v)      # Matching caps
        p_fixed[k] = np.real(v.conj().T @ result.lmp_masked[1] @ v)      # Fixed & tuning caps
        p_decouple[k] = np.real(v.conj().T @ result.lmp_masked[2] @ v)   # Decoupling cap

# Convert to percentage of forward power
power_data = np.column_stack([
    p_phantom, p_substrate, p_match, p_fixed, p_decouple, p_rad, p_ref
]) / forward_power * 100

# %% Create the plot
fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor('white')

# Create stacked area plot
ax.stackplot(phases, power_data.T, 
            labels=['Phantom', 'Substrate', 'C_Match', 'C_Tune', 'C_Decouple', 'Radiated', 'Reflected'],
            alpha=0.8)

ax.set_xlim(0, 180)
ax.set_ylim(0, 100)
ax.set_xlabel('Relative channel phase [°]')
ax.set_ylabel('Fraction of forward power [%]')
ax.legend(loc='lower left')
ax.set_title('Power balance of a 2-channel loop array under phase-only shimming')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% Print summary statistics
print("\n=== Power Balance Summary ===")
print(f"Forward power: {forward_power*1000:.1f} mW")
print(f"Number of phase points: {n_points}")
print(f"Phase range: {phases[0]:.0f}° to {phases[-1]:.0f}°")
print("\nPower fractions at 0° phase difference:")
print(f"  Reflected: {p_ref[0]/forward_power*100:.1f}%")
print(f"  Radiated: {p_rad[0]/forward_power*100:.1f}%")
print(f"  Phantom loss: {p_phantom[0]/forward_power*100:.1f}%")
print(f"  Substrate loss: {p_substrate[0]/forward_power*100:.1f}%")
print(f"  Matching caps: {p_match[0]/forward_power*100:.1f}%")
print(f"  Fixed/tuning caps: {p_fixed[0]/forward_power*100:.1f}%")
print(f"  Decoupling cap: {p_decouple[0]/forward_power*100:.1f}%")

print("\nPower fractions at 90° phase difference:")
idx_90 = n_points // 2
print(f"  Reflected: {p_ref[idx_90]/forward_power*100:.1f}%")
print(f"  Radiated: {p_rad[idx_90]/forward_power*100:.1f}%")
print(f"  Phantom loss: {p_phantom[idx_90]/forward_power*100:.1f}%")
print(f"  Substrate loss: {p_substrate[idx_90]/forward_power*100:.1f}%")
print(f"  Matching caps: {p_match[idx_90]/forward_power*100:.1f}%")
print(f"  Fixed/tuning caps: {p_fixed[idx_90]/forward_power*100:.1f}%")
print(f"  Decoupling cap: {p_decouple[idx_90]/forward_power*100:.1f}%")

# %% Power balance error analysis
print("\n=== Power Balance Error Analysis ===")

# Calculate total loss matrix (sum of all loss mechanisms)
# Masked matrices are subsets, so only use the full matrices
total_loss = (result.cpl + result.rad + result.mat + result.lmp)

# Power balance error matrix: Forward - Total_losses
power_error_matrix = result.fwd - total_loss

# Calculate eigenvalues to find worst-case error
eigenvalues = np.linalg.eigvals(power_error_matrix)

# Calculate error statistics
max_error_abs = np.max(np.abs(eigenvalues))
max_error_percent = (max_error_abs / forward_power) * 100

print(f"Power balance error matrix eigenvalues:")
for i, eig in enumerate(eigenvalues):
    print(f"  λ_{i+1}: {eig:.6f} W ({eig/forward_power*100:.3f}%)")

print(f"\nWorst-case power balance error: {max_error_abs:.6f} W ({max_error_percent:.3f}%)")

# Calculate relative error for each phase point
phase_errors = np.zeros(n_points)
for k in range(n_points):
    v = v_vectors[:, k]
    
    # Total power in = forward power for this excitation
    p_forward = np.real(v.conj().T @ result.fwd @ v)
    
    # Total power out = sum of all loss mechanisms (using full matrices only)
    p_coupling = np.real(v.conj().T @ result.cpl @ v)
    p_radiated = np.real(v.conj().T @ result.rad @ v)
    p_material = np.real(v.conj().T @ result.mat @ v)
    p_lumped = np.real(v.conj().T @ result.lmp @ v)
    
    p_total_loss = p_coupling + p_radiated + p_material + p_lumped
    
    # Power balance error
    phase_errors[k] = np.abs(p_forward - p_total_loss) / p_forward * 100

print(f"\nPhase-dependent power balance errors:")
print(f"  Mean error: {np.mean(phase_errors):.3f}%")
print(f"  Max error: {np.max(phase_errors):.3f}%")
print(f"  Min error: {np.min(phase_errors):.3f}%")
print(f"  Std error: {np.std(phase_errors):.3f}%")

# Verify that masked losses sum to total (as a sanity check)
if result.mat_masked:
    total_mat_from_masks = sum(result.mat_masked)
    mat_difference = np.max(np.abs(result.mat - total_mat_from_masks))
    print(f"\nMaterial loss partitioning check:")
    print(f"  Max difference between full matrix and sum of masks: {mat_difference:.2e}")

if result.lmp_masked:
    total_lmp_from_masks = sum(result.lmp_masked)
    lmp_difference = np.max(np.abs(result.lmp - total_lmp_from_masks))
    print(f"\nLumped loss partitioning check:")
    print(f"  Max difference between full matrix and sum of masks: {lmp_difference:.2e}")

# Check if error is acceptable (typically should be < 1%)
if max_error_percent < 1.0:
    print(f"\n✓ Power balance error ({max_error_percent:.3f}%) is within acceptable limits (<1%)")
else:
    print(f"\n⚠ Power balance error ({max_error_percent:.3f}%) exceeds typical limits (>1%)")
    print("  This may indicate issues with:")
    print("  - Grid resolution (too coarse)")
    print("  - Far-field offset (too small)")
    print("  - Boundary conditions")
    print("  - Field interpolation accuracy")
# %%
