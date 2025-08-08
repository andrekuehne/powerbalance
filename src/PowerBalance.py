"""
Power Balance Calculator for Multiport EM FDTD Simulations

This package calculates power balance matrices for electromagnetic FDTD simulations
on Yee grids, including forward power, coupling losses, material losses, lumped
element losses, and radiated power.
"""

from typing import Optional, List, Tuple, Union
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import RegularGridInterpolator

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
    """
    Calculate loss matrix from power radiated into far field.
    Matches MATLAB's calcRadiatedLossMatrix implementation.
    
    Args:
        e_field: Electric field array of shape [nPorts, nx, ny, nz, 3]
        h_field: Magnetic field array of shape [nPorts, nx, ny, nz, 3]
        x_line: x coordinate vector of Yee cell vertices (length nx+1)
        y_line: y coordinate vector of Yee cell vertices (length ny+1)
        z_line: z coordinate vector of Yee cell vertices (length nz+1)
        end_offset: Distance of integration box from outer boundaries in Yee cells
        
    Returns:
        Power correlation matrix of shape [nPorts, nPorts]
    """
    n_coils = e_field.shape[0]
    nx, ny, nz = e_field.shape[1:4]
    
    # Validate that end_offset is reasonable for the grid size
    if end_offset >= min(nx, ny, nz) // 2:
        raise ValueError(f"end_offset ({end_offset}) too large for grid size ({nx}x{ny}x{nz})")
    
    # Define grid lines and cell centers (matching MATLAB)
    dx = np.diff(x_line)
    dy = np.diff(y_line)
    dz = np.diff(z_line)
    
    cx_line = x_line[:-1] + dx / 2  # Cell centers in x
    cy_line = y_line[:-1] + dy / 2  # Cell centers in y
    cz_line = z_line[:-1] + dz / 2  # Cell centers in z
    
    # Storage for correlation matrices from each face
    matrices = []
    
    # Process each face direction and sign
    # X-directed faces
    matrices.extend(_process_x_faces(
        e_field, h_field, x_line, y_line, z_line,
        cx_line, cy_line, cz_line, dx, dy, dz,
        end_offset, n_coils
    ))
    
    # Y-directed faces
    matrices.extend(_process_y_faces(
        e_field, h_field, x_line, y_line, z_line,
        cx_line, cy_line, cz_line, dx, dy, dz,
        end_offset, n_coils
    ))
    
    # Z-directed faces
    matrices.extend(_process_z_faces(
        e_field, h_field, x_line, y_line, z_line,
        cx_line, cy_line, cz_line, dx, dy, dz,
        end_offset, n_coils
    ))
    
    # Sum all face contributions and enforce hermiticity
    if not matrices:
        return np.zeros((n_coils, n_coils), dtype=e_field.dtype)
    
    pcm = 0.5 * sum(matrices).T
    pcm = (pcm + pcm.conj().T) / 2
    
    return pcm


def _process_x_faces(e_field, h_field, x_line, y_line, z_line,
                     cx_line, cy_line, cz_line, dx, dy, dz,
                     end_offset, n_coils):
    """Process x-directed Yee faces (positive and negative)."""
    matrices = []
    
    # Surface area for x-faces
    sy, sz = np.meshgrid(
        dy[end_offset:len(dy)-end_offset],
        dz[end_offset:len(dz)-end_offset],
        indexing='ij'
    )
    surf_area_base = sy * sz
    
    # Positive x face (xp)
    x_face_idx = len(x_line) - 1 - end_offset
    if 0 <= x_face_idx < len(x_line):
        surf_area = surf_area_base.ravel()
        matrix = _calc_x_face_matrix(
            e_field, h_field, x_line, y_line, z_line,
            cx_line, cy_line, cz_line,
            x_face_idx, end_offset, n_coils, surf_area, 'positive'
        )
        if matrix is not None:
            matrices.append(matrix)
    
    # Negative x face (xm)
    x_face_idx = end_offset
    if 0 <= x_face_idx < len(x_line):
        surf_area = -surf_area_base.ravel()  # Negative for inward normal
        matrix = _calc_x_face_matrix(
            e_field, h_field, x_line, y_line, z_line,
            cx_line, cy_line, cz_line,
            x_face_idx, end_offset, n_coils, surf_area, 'negative'
        )
        if matrix is not None:
            matrices.append(matrix)
    
    return matrices


def _process_y_faces(e_field, h_field, x_line, y_line, z_line,
                     cx_line, cy_line, cz_line, dx, dy, dz,
                     end_offset, n_coils):
    """Process y-directed Yee faces (positive and negative)."""
    matrices = []
    
    # Surface area for y-faces
    sx, sz = np.meshgrid(
        dx[end_offset:len(dx)-end_offset],
        dz[end_offset:len(dz)-end_offset],
        indexing='ij'
    )
    surf_area_base = sx * sz
    
    # Positive y face (yp)
    y_face_idx = len(y_line) - 1 - end_offset
    if 0 <= y_face_idx < len(y_line):
        surf_area = surf_area_base.ravel()
        matrix = _calc_y_face_matrix(
            e_field, h_field, x_line, y_line, z_line,
            cx_line, cy_line, cz_line,
            y_face_idx, end_offset, n_coils, surf_area, 'positive'
        )
        if matrix is not None:
            matrices.append(matrix)
    
    # Negative y face (ym)
    y_face_idx = end_offset
    if 0 <= y_face_idx < len(y_line):
        surf_area = -surf_area_base.ravel()
        matrix = _calc_y_face_matrix(
            e_field, h_field, x_line, y_line, z_line,
            cx_line, cy_line, cz_line,
            y_face_idx, end_offset, n_coils, surf_area, 'negative'
        )
        if matrix is not None:
            matrices.append(matrix)
    
    return matrices


def _process_z_faces(e_field, h_field, x_line, y_line, z_line,
                     cx_line, cy_line, cz_line, dx, dy, dz,
                     end_offset, n_coils):
    """Process z-directed Yee faces (positive and negative)."""
    matrices = []
    
    # Surface area for z-faces
    sx, sy = np.meshgrid(
        dx[end_offset:len(dx)-end_offset],
        dy[end_offset:len(dy)-end_offset],
        indexing='ij'
    )
    surf_area_base = sx * sy
    
    # Positive z face (zp)
    z_face_idx = len(z_line) - 1 - end_offset
    if 0 <= z_face_idx < len(z_line):
        surf_area = surf_area_base.ravel()
        matrix = _calc_z_face_matrix(
            e_field, h_field, x_line, y_line, z_line,
            cx_line, cy_line, cz_line,
            z_face_idx, end_offset, n_coils, surf_area, 'positive'
        )
        if matrix is not None:
            matrices.append(matrix)
    
    # Negative z face (zm)
    z_face_idx = end_offset
    if 0 <= z_face_idx < len(z_line):
        surf_area = -surf_area_base.ravel()
        matrix = _calc_z_face_matrix(
            e_field, h_field, x_line, y_line, z_line,
            cx_line, cy_line, cz_line,
            z_face_idx, end_offset, n_coils, surf_area, 'negative'
        )
        if matrix is not None:
            matrices.append(matrix)
    
    return matrices


def _calc_x_face_matrix(e_field, h_field, x_line, y_line, z_line,
                        cx_line, cy_line, cz_line,
                        x_face_idx, end_offset, n_coils, surf_area, face_type):
    """Calculate correlation matrix for an x-directed face."""
    # Create evaluation points on the face
    x_face = x_line[x_face_idx]
    
    # Meshgrid for evaluation points
    mx, my, mz = np.meshgrid(
        x_face,
        cy_line[end_offset:len(cy_line)-end_offset],
        cz_line[end_offset:len(cz_line)-end_offset],
        indexing='ij'
    )
    
    # Flatten evaluation points
    eval_points = np.column_stack([
        mx.ravel(),
        my.ravel(),
        mz.ravel()
    ])
    
    # Storage for interpolated fields
    ey_interp = np.zeros((n_coils, len(eval_points)), dtype=e_field.dtype)
    ez_interp = np.zeros((n_coils, len(eval_points)), dtype=e_field.dtype)
    hy_interp = np.zeros((n_coils, len(eval_points)), dtype=h_field.dtype)
    hz_interp = np.zeros((n_coils, len(eval_points)), dtype=h_field.dtype)
    
    # Determine interpolation region based on face type
    if face_type == 'positive':
        x_idx_e = x_face_idx - 1  # E-field index
        x_range_e = slice(max(0, x_idx_e - 1), min(len(x_line) - 1, x_idx_e + 2))
        x_idx_h = x_face_idx - 1  # H-field center index
        x_range_h = slice(max(0, x_idx_h - 1), min(len(cx_line), x_idx_h + 2))
    else:  # negative
        x_idx_e = x_face_idx - 1
        x_range_e = slice(max(0, x_idx_e - 1), min(len(x_line) - 1, x_idx_e + 2))
        x_idx_h = x_face_idx - 1
        x_range_h = slice(max(0, x_idx_h - 1), min(len(cx_line), x_idx_h + 2))
    
    # Interpolate fields for each coil
    for coil in range(n_coils):
        # Ey interpolation (on x_line, cy_line, z_line[:-1])
        ey_interp[coil] = _interpolate_field_component(
            e_field[coil, x_range_e, :, :, 1],
            x_line[x_range_e], cy_line, z_line[:-1],
            eval_points
        )
        
        # Ez interpolation (on x_line, y_line[:-1], cz_line)
        ez_interp[coil] = _interpolate_field_component(
            e_field[coil, x_range_e, :, :, 2],
            x_line[x_range_e], y_line[:-1], cz_line,
            eval_points
        )
        
        # Hy interpolation (on cx_line, y_line[:-1], cz_line)
        hy_interp[coil] = _interpolate_field_component(
            h_field[coil, x_range_h, :, :, 1],
            cx_line[x_range_h], y_line[:-1], cz_line,
            eval_points
        )
        
        # Hz interpolation (on cx_line, cy_line, z_line[:-1])
        hz_interp[coil] = _interpolate_field_component(
            h_field[coil, x_range_h, :, :, 2],
            cx_line[x_range_h], cy_line, z_line[:-1],
            eval_points
        )
    
    # Calculate correlation matrix: Ey*surfArea @ Hz' - Ez*surfArea @ Hy'
    matrix = (
        (ey_interp * surf_area) @ hz_interp.conj().T -
        (ez_interp * surf_area) @ hy_interp.conj().T
    )
    
    return matrix


def _calc_y_face_matrix(e_field, h_field, x_line, y_line, z_line,
                        cx_line, cy_line, cz_line,
                        y_face_idx, end_offset, n_coils, surf_area, face_type):
    """Calculate correlation matrix for a y-directed face."""
    # Create evaluation points on the face
    y_face = y_line[y_face_idx]
    
    mx, my, mz = np.meshgrid(
        cx_line[end_offset:len(cx_line)-end_offset],
        y_face,
        cz_line[end_offset:len(cz_line)-end_offset],
        indexing='ij'
    )
    
    eval_points = np.column_stack([
        mx.ravel(),
        my.ravel(),
        mz.ravel()
    ])
    
    # Storage for interpolated fields
    ex_interp = np.zeros((n_coils, len(eval_points)), dtype=e_field.dtype)
    ez_interp = np.zeros((n_coils, len(eval_points)), dtype=e_field.dtype)
    hx_interp = np.zeros((n_coils, len(eval_points)), dtype=h_field.dtype)
    hz_interp = np.zeros((n_coils, len(eval_points)), dtype=h_field.dtype)
    
    # Determine interpolation region
    if face_type == 'positive':
        y_idx_e = y_face_idx - 1
        y_range_e = slice(max(0, y_idx_e - 1), min(len(y_line) - 1, y_idx_e + 2))
        y_idx_h = y_face_idx - 1
        y_range_h = slice(max(0, y_idx_h - 1), min(len(cy_line), y_idx_h + 2))
    else:
        y_idx_e = y_face_idx - 1
        y_range_e = slice(max(0, y_idx_e - 1), min(len(y_line) - 1, y_idx_e + 2))
        y_idx_h = y_face_idx - 1
        y_range_h = slice(max(0, y_idx_h - 1), min(len(cy_line), y_idx_h + 2))
    
    # Interpolate fields for each coil
    for coil in range(n_coils):
        # Ex interpolation (on cx_line, y_line, z_line[:-1])
        ex_interp[coil] = _interpolate_field_component(
            e_field[coil, :, y_range_e, :, 0],
            cx_line, y_line[y_range_e], z_line[:-1],
            eval_points
        )
        
        # Ez interpolation (on x_line[:-1], y_line, cz_line)
        ez_interp[coil] = _interpolate_field_component(
            e_field[coil, :, y_range_e, :, 2],
            x_line[:-1], y_line[y_range_e], cz_line,
            eval_points
        )
        
        # Hx interpolation (on x_line[:-1], cy_line, cz_line)
        hx_interp[coil] = _interpolate_field_component(
            h_field[coil, :, y_range_h, :, 0],
            x_line[:-1], cy_line[y_range_h], cz_line,
            eval_points
        )
        
        # Hz interpolation (on cx_line, cy_line, z_line[:-1])
        hz_interp[coil] = _interpolate_field_component(
            h_field[coil, :, y_range_h, :, 2],
            cx_line, cy_line[y_range_h], z_line[:-1],
            eval_points
        )
    
    # Calculate correlation matrix: Ez*surfArea @ Hx' - Ex*surfArea @ Hz'
    matrix = (
        (ez_interp * surf_area) @ hx_interp.conj().T -
        (ex_interp * surf_area) @ hz_interp.conj().T
    )
    
    return matrix


def _calc_z_face_matrix(e_field, h_field, x_line, y_line, z_line,
                        cx_line, cy_line, cz_line,
                        z_face_idx, end_offset, n_coils, surf_area, face_type):
    """Calculate correlation matrix for a z-directed face."""
    # Create evaluation points on the face
    z_face = z_line[z_face_idx]
    
    mx, my, mz = np.meshgrid(
        cx_line[end_offset:len(cx_line)-end_offset],
        cy_line[end_offset:len(cy_line)-end_offset],
        z_face,
        indexing='ij'
    )
    
    eval_points = np.column_stack([
        mx.ravel(),
        my.ravel(),
        mz.ravel()
    ])
    
    # Storage for interpolated fields
    ex_interp = np.zeros((n_coils, len(eval_points)), dtype=e_field.dtype)
    ey_interp = np.zeros((n_coils, len(eval_points)), dtype=e_field.dtype)
    hx_interp = np.zeros((n_coils, len(eval_points)), dtype=h_field.dtype)
    hy_interp = np.zeros((n_coils, len(eval_points)), dtype=h_field.dtype)
    
    # Determine interpolation region
    if face_type == 'positive':
        z_idx_e = z_face_idx - 1
        z_range_e = slice(max(0, z_idx_e - 1), min(len(z_line) - 1, z_idx_e + 2))
        z_idx_h = z_face_idx - 1
        z_range_h = slice(max(0, z_idx_h - 1), min(len(cz_line), z_idx_h + 2))
    else:
        z_idx_e = z_face_idx - 1
        z_range_e = slice(max(0, z_idx_e - 1), min(len(z_line) - 1, z_idx_e + 2))
        z_idx_h = z_face_idx - 1
        z_range_h = slice(max(0, z_idx_h - 1), min(len(cz_line), z_idx_h + 2))
    
    # Interpolate fields for each coil
    for coil in range(n_coils):
        # Ex interpolation (on cx_line, y_line[:-1], z_line)
        ex_interp[coil] = _interpolate_field_component(
            e_field[coil, :, :, z_range_e, 0],
            cx_line, y_line[:-1], z_line[z_range_e],
            eval_points
        )
        
        # Ey interpolation (on x_line[:-1], cy_line, z_line)
        ey_interp[coil] = _interpolate_field_component(
            e_field[coil, :, :, z_range_e, 1],
            x_line[:-1], cy_line, z_line[z_range_e],
            eval_points
        )
        
        # Hx interpolation (on x_line[:-1], cy_line, cz_line)
        hx_interp[coil] = _interpolate_field_component(
            h_field[coil, :, :, z_range_h, 0],
            x_line[:-1], cy_line, cz_line[z_range_h],
            eval_points
        )
        
        # Hy interpolation (on cx_line, y_line[:-1], cz_line)
        hy_interp[coil] = _interpolate_field_component(
            h_field[coil, :, :, z_range_h, 1],
            cx_line, y_line[:-1], cz_line[z_range_h],
            eval_points
        )
    
    # Calculate correlation matrix: Ex*surfArea @ Hy' - Ey*surfArea @ Hx'
    matrix = (
        (ex_interp * surf_area) @ hy_interp.conj().T -
        (ey_interp * surf_area) @ hx_interp.conj().T
    )
    
    return matrix


def _interpolate_field_component(field_data, x_grid, y_grid, z_grid, eval_points):
    """
    Interpolate a field component to evaluation points.
    Uses scipy's RegularGridInterpolator to match MATLAB's griddedInterpolant.
    
    Args:
        field_data: 3D field data array
        x_grid, y_grid, z_grid: Grid coordinate vectors
        eval_points: Points where field should be evaluated
        
    Returns:
        Interpolated field values at evaluation points
    """
    # Handle edge cases where grid is too small
    if len(x_grid) < 2 or len(y_grid) < 2 or len(z_grid) < 2:
        return np.zeros(len(eval_points), dtype=field_data.dtype)
    
    # Ensure field data matches grid dimensions
    expected_shape = (len(x_grid), len(y_grid), len(z_grid))
    if field_data.shape != expected_shape:
        # Reshape or pad as needed
        reshaped = np.zeros(expected_shape, dtype=field_data.dtype)
        min_shape = tuple(min(s1, s2) for s1, s2 in zip(field_data.shape, expected_shape))
        reshaped[:min_shape[0], :min_shape[1], :min_shape[2]] = \
            field_data[:min_shape[0], :min_shape[1], :min_shape[2]]
        field_data = reshaped
    
    # Create interpolator
    # Use 'linear' method to match MATLAB's default behavior
    # bounds_error=False and fill_value=0 to handle edge cases
    interpolator = RegularGridInterpolator(
        (x_grid, y_grid, z_grid),
        field_data,
        method='linear',
        bounds_error=False,
        fill_value=0.0
    )
    
    # Evaluate at points
    return interpolator(eval_points)

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