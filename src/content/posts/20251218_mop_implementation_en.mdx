---
title: 【Computational Chemistry】 New Features in MultiOptPy v1.20.3
published: 2025-12-18
tags: [MultiOptPy, python]
category: Computational Chemistry
draft: false
---
Last Updated: 2025-12-18

### 1. Mode-Following RS-I-RFO (MF-RS-I-RFO)

#### Overview
This is an algorithm for Saddle Point Optimization that proceeds by tracking a specific Hessian mode. Based on the conventional RS-I-RFO method, it incorporates tracking logic to prevent mode swapping between steps.

#### Implementation Details
*   **Mode Tracking Strategy:** Uses the `ModeFollowing` class to calculate the overlap with the Hessian eigenvectors of the previous step.
*   **Mass-Weighted Overlap (MWO):** If an atom list is provided, the projection is performed in the mass-weighted coordinate system rather than Cartesian coordinates, achieving physically valid tracking.
*   **Adaptive Update (EMA):** Uses Exponential Moving Average (EMA) to dynamically update the reference vector, allowing the method to follow mode rotation.
*   **Gradient Bias:** Adds the overlap in the direction of the current force (gradient) to the score, prioritizing the selection of energetically significant modes.

#### Command Line Arguments (`-opt`)
Detailed parameters are specified using colons `:` within the `-opt` argument string.

*   **Basic Format:** `-opt mf_rsirfo:<target_index>:<ema_val>:<grad_val>`
*   **Parameters:**
    *   `target_index`: The index of the eigenvalue to track (starting from 0).
    *   `ema<val>`: Update rate (0.0 to 1.0). 1.0 for full replacement (adaptive), 0.0 for fixed.
    *   `grad<val>`: Weighting for the gradient direction.

#### Usage Example:
```bash
# Execute transition state geometry optimization tracking mass-weighted mode 1, 
# with an adaptive update rate of 0.5 and a gradient bias of 0.3.
python optmain.py input.xyz -opt mwmf_rsirfo_fsb:1:ema0.5:grad0.3 -fc 5 -order 1 -freq -tcc
```

### 2. Constrained RS-I-RFO (C-RS-I-RFO)

#### Overview
A method for performing saddle point searches using RS-I-RFO while imposing geometric constraints (bond distances, angles, etc.).

#### Implementation Details
*   **Subspace Projection:** Uses Singular Value Decomposition (SVD) to construct a "null space basis" orthogonal to the imposed constraints.
*   **Constrained Optimization:** Projects the full-space gradient and Hessian onto this subspace and calculates the RFO step within it. This calculates the movement vector to the saddle point while maintaining the constraints.
*   **SHAKE-like Correction:** Includes geometric adjustment logic to correct deviations in constraints due to numerical errors.

#### Command Line Arguments (`-opt`, `-pc`)
Specify `crsirfo` family methods with `-opt` and define constraints with `-pc`.

#### Usage Example:
```bash
# Geometry optimization while fixing the distance between atoms 1 and 2, and the angle of atoms 2-3-4.
python optmain.py input.xyz -opt crsirfo_block_fsb -pc bond 1,2 bend 2,3,4
```

### 3. Integration of Multiple PES Information and Model Function Optimization (BITSS, etc.)

#### Overview
A feature to perform geometry optimization on an "effective potential" constructed by combining energies or gradients from multiple electronic states (PES).
Specifically, version 1.20.3 implements the Binary-Image Transition State Search (BITSS) method.

#### Implementation Details
*   **Independent Calculator Instances:** `ModelFunctionHandler` creates independent directories and calculators for State1 and State2 to prevent interference between states.
*   **BITSS Implementation:**
    *   **6N-Dimensional Expansion:** Concatenates two structures (Image 1 & 2) and treats them as a $6N \times 6N$ Hessian and a $2N$ atom system.
    *   **Second Derivative (Hessian):** Mathematically derives the second derivatives for distance and energy equality constraint terms and implements them as the `calc_hess` method.

#### Command Line Arguments (`-mf`)
Use the `-mf` argument to specify the type of model function and accompanying parameters (the file path of the reference structure in the case of BITSS).

#### Usage Example:
```bash
# Search for a pathway connecting two points using the BITSS method, with target.xyz as the target structure.
python optmain.py start.xyz -mf bitss target.xyz
```

```bash
# Minimum energy seam of crossing (MESX) using Seam Model Function (between charge 0, multiplicity 1 and 3 states)
python optmain.py input.xyz -mf opt_meci 0 3 -elec 0 -spin 3
```

#### Note: Internal Data Structure Changes
*   **OptimizationState:** In BITSS mode, `element_list` and `geometry` are automatically expanded to double size (2N).
*   **Hessian Integration:** `Model_hess` (derived from PES) and `bias_hessian` (derived from constraint/penalty terms) are managed separately and summed just before being passed to the optimizer.


#### References
- _J. Chem. Phys._ 157, 124107 (2022)
- _J. Am. Chem. Soc._ 2015, 137, 10, 3433–3445
