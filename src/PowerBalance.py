"""
Power Balance Calculator for Multiport EM FDTD Simulations

This package calculates power balance matrices for electromagnetic FDTD simulations
on Yee grids, including forward power, coupling losses, material losses, lumped
element losses, and radiated power.
"""

from typing import Optional, List, Tuple, Union
from dataclasses import dataclass
import numpy as np


@dataclass
class PowerBalanceResult:
    """Container for power balance calculation results."""
    fwd: np.ndarray  # Forward power matrix
    cpl: np.ndarray  # Coupling loss matrix
    lmp: np.ndarray  # Lumped element loss matrix
    mat: np.ndarray  # Material loss matrix
    rad: np.ndarray  # Radiated power matrix
    lmp_masked: Optional[List[np.ndarray]] = None  # Masked lumped losses
    mat_masked: Optional[List[np.ndarray]] = None  # Masked material losses


class PowerBalanceCalculator:
    """
    Calculator for power balance matrices in multiport EM FDTD simulations.
    
    This class orchestrates the calculation of various power loss mechanisms
    in electromagnetic simulations on Yee grids.
    """
    
    def __init__(self, far_field_offset: int = 4):
        """
        Initialize calculator.
        
        Args:
            far_field_offset: Distance from outer field boundaries for 
                            Poynting flux integration
        """
        self.far_field_offset = far_field_offset
    
    def calculate_power_balance(
        self,
        e_field: np.ndarray,
        h_field: np.ndarray,
        sigma: np.ndarray,
        grid_vectors: Tuple[np.ndarray, np.ndarray, np.ndarray],
        u_matrix: np.ndarray,
        i_matrix: np.ndarray,
        s_matrix: np.ndarray,
        forward_power: float,
        sigma_masks: Optional[List[np.ndarray]] = None,
        lumped_masks: Optional[List[np.ndarray]] = None
    ) -> PowerBalanceResult:
        """
        Calculate power balance matrices for multiport EM FDTD simulation.
        
        Args:
            e_field: Electric field array of shape [nPorts, nx, ny, nz, 3]
            h_field: Magnetic field array of shape [nPorts, nx, ny, nz, 3]
            sigma: Electric conductivity array of shape [nx, ny, nz, 3]
            grid_vectors: Tuple of (x_line, y_line, z_line) coordinate arrays
            u_matrix: Lumped element voltages [nLumpedElements, nPorts]
            i_matrix: Lumped element currents [nLumpedElements, nPorts]
            s_matrix: Scattering matrix [nPorts, nPorts]
            forward_power: Forward power scalar
            sigma_masks: Optional list of conductivity masks
            lumped_masks: Optional list of lumped element masks
            
        Returns:
            PowerBalanceResult containing all calculated matrices
        """
        # Validate inputs
        self._validate_inputs(e_field, h_field, sigma, grid_vectors, 
                            u_matrix, i_matrix, s_matrix, forward_power,
                            sigma_masks, lumped_masks)
        
        n_ports = e_field.shape[0]
        x_line, y_line, z_line = grid_vectors
        
        # Calculate forward power matrix
        fwd = calc_forward_power_matrix(n_ports, forward_power)
        
        # Calculate coupling loss matrix
        cpl = calc_coupling_loss_matrix(s_matrix, forward_power)
        
        # Calculate lumped element loss matrix
        lmp = calc_lumped_loss_matrix(u_matrix, i_matrix)
        
        # Calculate masked lumped losses if masks provided
        lmp_masked = None
        if lumped_masks is not None:
            lmp_masked = []
            for mask in lumped_masks:
                curr_u = u_matrix * mask[:, np.newaxis]
                curr_i = i_matrix * mask[:, np.newaxis]
                lmp_masked.append(calc_lumped_loss_matrix(curr_u, curr_i))
        
        # Calculate material loss matrix
        mat = calc_material_loss_matrix(e_field, sigma, x_line, y_line, z_line)
        
        # Calculate masked material losses if masks provided
        mat_masked = None
        if sigma_masks is not None:
            mat_masked = []
            for mask in sigma_masks:
                masked_sigma = sigma * mask
                mat_masked.append(calc_material_loss_matrix(
                    e_field, masked_sigma, x_line, y_line, z_line))
        
        # Calculate radiated loss matrix
        rad = calc_radiated_loss_matrix(e_field, h_field, x_line, y_line, 
                                      z_line, self.far_field_offset)
        
        return PowerBalanceResult(
            fwd=fwd, cpl=cpl, lmp=lmp, mat=mat, rad=rad,
            lmp_masked=lmp_masked, mat_masked=mat_masked
        )
    
    def _validate_inputs(self, e_field, h_field, sigma, grid_vectors,
                        u_matrix, i_matrix, s_matrix, forward_power,
                        sigma_masks, lumped_masks):
        """Validate all input parameters."""
        # E-field validation
        if e_field.ndim != 5:
            raise ValueError('E-field must have 5 dimensions [nPorts, nx, ny, nz, 3]')
        if e_field.shape[-1] != 3:
            raise ValueError('E-field last dimension must be 3')
        if not (np.issubdtype(e_field.dtype, np.floating) or np.issubdtype(e_field.dtype, np.complexfloating)):
            raise TypeError('E-field must be floating point or complex floating point type')
        
        sz_fields = e_field.shape
        grid_dims = sz_fields[1:4]
        n_ports = sz_fields[0]
        
        # H-field validation
        if h_field.shape != e_field.shape:
            raise ValueError('H-field must have same shape as E-field')
        if not (np.issubdtype(h_field.dtype, np.floating) or np.issubdtype(h_field.dtype, np.complexfloating)):
            raise TypeError('H-field must be floating point or complex floating point type')
        
        # Sigma validation
        if sigma.ndim != 4:
            raise ValueError('Sigma must have 4 dimensions [nx, ny, nz, 3]')
        if sigma.shape != (*grid_dims, 3):
            raise ValueError('Sigma size does not match grid dimensions')
        if not np.issubdtype(sigma.dtype, np.floating):
            raise TypeError('Sigma must be floating point type')
        
        # Grid vectors validation
        if len(grid_vectors) != 3:
            raise ValueError('Three grid vectors required')
        for i, vec in enumerate(grid_vectors):
            if not isinstance(vec, np.ndarray) or vec.ndim != 1:
                raise ValueError(f'Grid vector {i} must be 1D array')
            if not np.issubdtype(vec.dtype, np.floating):
                raise TypeError(f'Grid vector {i} must be floating point type')
        
        sz_grid = [len(vec) for vec in grid_vectors]
        if tuple(grid_dims) != tuple(np.array(sz_grid) - 1):
            raise ValueError('Grid size incompatible with field size')
        
        # Matrix validations
        if u_matrix.ndim != 2:
            raise ValueError('Voltage matrix must be 2D')
        if u_matrix.shape[1] != n_ports:
            raise ValueError('Voltage matrix columns must match number of ports')
        if not (np.issubdtype(u_matrix.dtype, np.floating) or np.issubdtype(u_matrix.dtype, np.complexfloating)):
            raise TypeError('Voltage matrix must be floating point or complex floating point type')
        
        if i_matrix.shape != u_matrix.shape:
            raise ValueError('Current matrix must have same shape as voltage matrix')
        if not (np.issubdtype(i_matrix.dtype, np.floating) or np.issubdtype(i_matrix.dtype, np.complexfloating)):
            raise TypeError('Current matrix must be floating point or complex floating point type')
        
        if s_matrix.ndim != 2 or s_matrix.shape[0] != s_matrix.shape[1]:
            raise ValueError('S-matrix must be square')
        if s_matrix.shape[0] != n_ports:
            raise ValueError('S-matrix size must match number of ports')
        if not (np.issubdtype(s_matrix.dtype, np.floating) or np.issubdtype(s_matrix.dtype, np.complexfloating)):
            raise TypeError('S-matrix must be floating point or complex floating point type')
        
        # Forward power validation
        if not np.isscalar(forward_power) or not np.isreal(forward_power):
            raise ValueError('Forward power must be real scalar')
        
        # Mask validations
        if sigma_masks is not None:
            for i, mask in enumerate(sigma_masks):
                if mask.shape != sigma.shape:
                    raise ValueError(f'Sigma mask {i} shape mismatch')
                if mask.dtype != bool:
                    raise TypeError(f'Sigma mask {i} must be boolean')
        
        if lumped_masks is not None:
            n_lumped = u_matrix.shape[0]
            for i, mask in enumerate(lumped_masks):
                if mask.shape != (n_lumped,):
                    raise ValueError(f'Lumped mask {i} shape mismatch')
                if mask.dtype != bool:
                    raise TypeError(f'Lumped mask {i} must be boolean')


def calc_forward_power_matrix(n_ports: int, forward_power: float) -> np.ndarray:
    """Calculate forward power matrix from forward power and number of ports."""
    return np.eye(n_ports) * forward_power


def calc_coupling_loss_matrix(s_matrix: np.ndarray, forward_power: float) -> np.ndarray:
    """Calculate coupling loss matrix from scattering matrix and forward power."""
    # For complex S-parameters: |S|² = S* @ S
    return s_matrix.conj().T @ s_matrix * forward_power


def calc_lumped_loss_matrix(u_matrix: np.ndarray, i_matrix: np.ndarray) -> np.ndarray:
    """Calculate loss matrix from lumped element voltages and currents."""
    # For complex phasors, power = Re(V * I*) = Re((I* * V + V* * I)/2)
    pcm = (i_matrix.conj().T @ u_matrix + u_matrix.conj().T @ i_matrix) / 4
    pcm = (pcm + pcm.conj().T) / 2  # Enforce hermiticity
    return pcm


def calc_material_loss_matrix(
    e_field: np.ndarray,
    sigma: np.ndarray,
    x_line: np.ndarray,
    y_line: np.ndarray,
    z_line: np.ndarray
) -> np.ndarray:
    """Calculate loss matrix from electric fields in lossy materials."""
    # Generate grid variables
    dx = np.diff(x_line)
    dy = np.diff(y_line)
    dz = np.diff(z_line)
    
    cx_line = x_line[:-1] + dx / 2
    cy_line = y_line[:-1] + dy / 2
    cz_line = z_line[:-1] + dz / 2
    
    dcx = np.diff(cx_line)
    dcy = np.diff(cy_line)
    dcz = np.diff(cz_line)
    
    # Generate grids for integration deltas
    mxx, myx, mzx = np.meshgrid(dx, dcy, dcz, indexing='ij')
    vx = mxx * myx * mzx
    
    mxy, myy, mzy = np.meshgrid(dcx, dy, dcz, indexing='ij')
    vy = mxy * myy * mzy
    
    mxz, myz, mzz = np.meshgrid(dcx, dcy, dz, indexing='ij')
    vz = mxz * myz * mzz
    
    # Pad grid with zeros for unusable boundary regions
    vx = np.pad(vx, ((0, 0), (1, 0), (1, 0)), constant_values=0)
    vy = np.pad(vy, ((1, 0), (0, 0), (1, 0)), constant_values=0)
    vz = np.pad(vz, ((1, 0), (1, 0), (0, 0)), constant_values=0)
    
    # Concatenate to full 3D vector grid
    v = np.stack([vx, vy, vz], axis=3)
    
    # Find points with nonzero conductivity
    valid_edges = sigma > 0
    
    # Only use cell edges with conductivity
    e_valid = e_field[:, valid_edges]  # Shape: [nPorts, nValidEdges]
    s_valid = sigma[valid_edges]  # Shape: [nValidEdges]
    v_valid = v[valid_edges]  # Shape: [nValidEdges]
    
    # Calculate power correlation matrix
    # For complex fields: P = 0.5 * Re(E * conj(E) * sigma * volume)
    weighted_e = e_valid * (s_valid * v_valid)[np.newaxis, :]
    pcm = 0.5 * (weighted_e @ e_valid.conj().T).T
    pcm = (pcm + pcm.conj().T) / 2  # Enforce hermiticity
    
    return pcm


def calc_radiated_loss_matrix(
    e_field: np.ndarray,
    h_field: np.ndarray,
    x_line: np.ndarray,
    y_line: np.ndarray,
    z_line: np.ndarray,
    end_offset: int
) -> np.ndarray:
    """Calculate loss matrix from power radiated into far field."""
    n_coils = e_field.shape[0]
    nx, ny, nz = e_field.shape[1:4]
    
    # Validate that end_offset is reasonable for the grid size
    if end_offset >= min(nx, ny, nz) // 2:
        raise ValueError(f"end_offset ({end_offset}) too large for grid size ({nx}x{ny}x{nz})")
    
    # Define gridlines
    dx = np.diff(x_line)
    dy = np.diff(y_line)
    dz = np.diff(z_line)
    
    # Calculate Poynting vector on boundary surfaces
    surfaces = _calc_poynting_surfaces_simple(
        e_field, h_field, dx, dy, dz, end_offset, n_coils
    )
    
    # Calculate power correlation matrices for each surface - following MATLAB pattern
    matrices = []
    for direction in ['x', 'y', 'z']:
        for sign in ['p', 'm']:
            surf_key = f"{direction}{sign}"
            if surf_key not in surfaces:
                continue
                
            surf = surfaces[surf_key]
            surf_area = surf['surfArea'].ravel()
            
            # Extract and reshape field components 
            if direction == 'x':
                # For x-face: Poynting = Ey * Hz* - Ez * Hy*
                ey_flat = surf['ey'].reshape(n_coils, -1)
                ez_flat = surf['ez'].reshape(n_coils, -1)
                hy_flat = surf['hy'].reshape(n_coils, -1)
                hz_flat = surf['hz'].reshape(n_coils, -1)
                
                # Calculate correlation matrix: bsxfun(@times, Ey, surfArea) * Hz' - bsxfun(@times, Ez, surfArea) * Hy'
                matrix = ((ey_flat * surf_area) @ hz_flat.conj().T - 
                         (ez_flat * surf_area) @ hy_flat.conj().T)
                
            elif direction == 'y':
                # For y-face: Poynting = Ez * Hx* - Ex * Hz*
                ex_flat = surf['ex'].reshape(n_coils, -1)
                ez_flat = surf['ez'].reshape(n_coils, -1)
                hx_flat = surf['hx'].reshape(n_coils, -1)
                hz_flat = surf['hz'].reshape(n_coils, -1)
                
                matrix = ((ez_flat * surf_area) @ hx_flat.conj().T - 
                         (ex_flat * surf_area) @ hz_flat.conj().T)
                
            else:  # direction == 'z'
                # For z-face: Poynting = Ex * Hy* - Ey * Hx*
                ex_flat = surf['ex'].reshape(n_coils, -1)
                ey_flat = surf['ey'].reshape(n_coils, -1)
                hx_flat = surf['hx'].reshape(n_coils, -1)
                hy_flat = surf['hy'].reshape(n_coils, -1)
                
                matrix = ((ex_flat * surf_area) @ hy_flat.conj().T - 
                         (ey_flat * surf_area) @ hx_flat.conj().T)
            
            matrices.append(matrix)
    
    # Sum all surface contributions
    if not matrices:
        return np.zeros((n_coils, n_coils), dtype=complex)
    
    # The factor of 0.5 and transpose match the MATLAB implementation
    pcm = 0.5 * sum(matrices).T
    pcm = (pcm + pcm.conj().T) / 2  # Enforce hermiticity for complex result
    
    return pcm  # Return real part for power matrix


def _calc_poynting_surfaces_simple(e_field, h_field, dx, dy, dz, end_offset, n_coils):
    """Calculate Poynting vector on boundary surfaces using simple H-field averaging."""
    surfaces = {}
    nx, ny, nz = e_field.shape[1:4]
    
    # Helper function for safe averaging
    def safe_avg(field, idx1, idx2, axis):
        """Safely average two field slices along an axis."""
        if idx1 < 0 or idx2 >= field.shape[axis + 1]:  # +1 because first dim is n_coils
            return field.take(max(0, min(idx1, idx2)), axis=axis + 1)
        
        f1 = field.take(idx1, axis=axis + 1)
        f2 = field.take(idx2, axis=axis + 1)
        return (f1 + f2) / 2
    
    # X-directed faces
    for sign, x_idx in [('p', nx - 1 - end_offset), ('m', end_offset)]:
        if x_idx < 0 or x_idx >= nx:
            continue
            
        y_start, y_end = end_offset, ny - end_offset
        z_start, z_end = end_offset, nz - end_offset
        
        if y_start >= y_end or z_start >= z_end:
            continue
        
        # Surface area
        sy, sz = np.meshgrid(dy[y_start:y_end], dz[z_start:z_end], indexing='ij')
        surf_area = sy * sz
        if sign == 'm':
            surf_area = -surf_area
        
        # E field components (already at correct locations)
        ey = e_field[:, x_idx, y_start:y_end, z_start:z_end, 1]
        ez = e_field[:, x_idx, y_start:y_end, z_start:z_end, 2]
        
        # H field components (need interpolation from cell centers to face centers)
        # For x-face: need Hy and Hz interpolated to x_line locations
        hy = safe_avg(h_field[:, :, y_start:y_end, z_start:z_end, 1], x_idx-1, x_idx, 0)
        hz = safe_avg(h_field[:, :, y_start:y_end, z_start:z_end, 2], x_idx-1, x_idx, 0)
        
        surfaces[f'x{sign}'] = {
            'surfArea': surf_area,
            'ey': ey,
            'ez': ez,
            'hy': hy,
            'hz': hz,
        }
    
    # Y-directed faces
    for sign, y_idx in [('p', ny - 1 - end_offset), ('m', end_offset)]:
        if y_idx < 0 or y_idx >= ny:
            continue
            
        x_start, x_end = end_offset, nx - end_offset
        z_start, z_end = end_offset, nz - end_offset
        
        if x_start >= x_end or z_start >= z_end:
            continue
        
        sx, sz = np.meshgrid(dx[x_start:x_end], dz[z_start:z_end], indexing='ij')
        surf_area = sx * sz
        if sign == 'm':
            surf_area = -surf_area
        
        # E field components
        ex = e_field[:, x_start:x_end, y_idx, z_start:z_end, 0]
        ez = e_field[:, x_start:x_end, y_idx, z_start:z_end, 2]
        
        # H field components (interpolated to y-face)
        hx = safe_avg(h_field[:, x_start:x_end, :, z_start:z_end, 0], y_idx-1, y_idx, 1)
        hz = safe_avg(h_field[:, x_start:x_end, :, z_start:z_end, 2], y_idx-1, y_idx, 1)
        
        surfaces[f'y{sign}'] = {
            'surfArea': surf_area,
            'ex': ex,
            'ez': ez,
            'hx': hx,
            'hz': hz,
        }
    
    # Z-directed faces
    for sign, z_idx in [('p', nz - 1 - end_offset), ('m', end_offset)]:
        if z_idx < 0 or z_idx >= nz:
            continue
            
        x_start, x_end = end_offset, nx - end_offset
        y_start, y_end = end_offset, ny - end_offset
        
        if x_start >= x_end or y_start >= y_end:
            continue
        
        sx, sy = np.meshgrid(dx[x_start:x_end], dy[y_start:y_end], indexing='ij')
        surf_area = sx * sy
        if sign == 'm':
            surf_area = -surf_area
        
        # E field components
        ex = e_field[:, x_start:x_end, y_start:y_end, z_idx, 0]
        ey = e_field[:, x_start:x_end, y_start:y_end, z_idx, 1]
        
        # H field components (interpolated to z-face)
        hx = safe_avg(h_field[:, x_start:x_end, y_start:y_end, :, 0], z_idx-1, z_idx, 2)
        hy = safe_avg(h_field[:, x_start:x_end, y_start:y_end, :, 1], z_idx-1, z_idx, 2)
        
        surfaces[f'z{sign}'] = {
            'surfArea': surf_area,
            'ex': ex,
            'ey': ey,
            'hx': hx,
            'hy': hy,
        }
    
    return surfaces


# Example usage
if __name__ == "__main__":
    # Example usage with dummy data
    n_ports = 2
    nx, ny, nz = 20, 20, 20
    
    # Create dummy field data (complex for time-harmonic phasors)
    e_field = (np.random.randn(n_ports, nx, ny, nz, 3) + 
               1j * np.random.randn(n_ports, nx, ny, nz, 3)).astype(np.complex64)
    h_field = (np.random.randn(n_ports, nx, ny, nz, 3) + 
               1j * np.random.randn(n_ports, nx, ny, nz, 3)).astype(np.complex64)
    sigma = np.random.rand(nx, ny, nz, 3).astype(np.float32) * 0.1
    
    # Create grid vectors
    x_line = np.linspace(0, 1, nx + 1)
    y_line = np.linspace(0, 1, ny + 1)
    z_line = np.linspace(0, 1, nz + 1)
    grid_vectors = (x_line, y_line, z_line)
    
    # Create dummy lumped element data (complex)
    n_lumped = 5
    u_matrix = (np.random.randn(n_lumped, n_ports) + 
                1j * np.random.randn(n_lumped, n_ports)).astype(np.complex64)
    i_matrix = (np.random.randn(n_lumped, n_ports) + 
                1j * np.random.randn(n_lumped, n_ports)).astype(np.complex64)
    
    # Create S-matrix (complex)
    s_matrix = (np.random.randn(n_ports, n_ports) + 
                1j * np.random.randn(n_ports, n_ports)).astype(np.complex64) * 0.1
    forward_power = 1.0
    
    # Calculate power balance
    calculator = PowerBalanceCalculator(far_field_offset=2)  # Smaller offset for test
    result = calculator.calculate_power_balance(
        e_field, h_field, sigma, grid_vectors,
        u_matrix, i_matrix, s_matrix, forward_power
    )
    
    print("Power Balance Calculation Complete!")
    print(f"Forward power matrix shape: {result.fwd.shape}")
    print(f"Coupling loss matrix shape: {result.cpl.shape}")
    print(f"Lumped loss matrix shape: {result.lmp.shape}")
    print(f"Material loss matrix shape: {result.mat.shape}")
    print(f"Radiated power matrix shape: {result.rad.shape}")