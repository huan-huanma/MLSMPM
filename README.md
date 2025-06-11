# MLSMPM
MLSMPM for large deformation problems
# 1D Elastic Bar Vibration Simulation using MLS-MPM at MLSMPM/master/
This MATLAB script implements the **Moving Least Squares Material Point Method (MLS-MPM)** for simulating the transient vibration of a **1D elastic bar**. The approach combines meshfree MLS shape functions with MPM.

---

##  Features

- Simulates a **1D bar** undergoing elastic vibration.
- Uses MLS-based shape functions (`quartic_spline`) for **improved accuracy** over standard MPM.
- Computes and compares:
  - **Center-of-mass velocity** history.
  - **Kinetic, strain, and total energy** evolution.
  - **L2 error** between numerical and analytical solutions.
- Visualizes:
  - Velocity history (numerical vs analytical).
  - Absolute error.
  - Energy evolution.

---

##  Parameters

- **Material properties:**
  - Young's modulus \( E=100 \)
  - Density \( \rho=1 \)
  - Initial velocity amplitude \( v_0=0.1 \)
- **Simulation domain:** \( L=1 \)
- **MLS parameters:**
  - Shape: quartic spline weight function.
  - Support domain radius scaling \( d_\text{max}=2.5 \).
- **Time step:** adaptive based on CFL condition.
- **Mode number:** 1 (first mode vibration).

---

##  How to Run

1. Clone/download the repository.
2. Ensure the following helper functions are available:
   - `buildGrid1D.m`: Build the 1D computational grid.
   - `buildParticlesGeom.m`: Initialize particles.
   - `defineSupport.m`: Determine MLS support domain.
   - `mlsQuadricBasis1D.m`: Compute 1D MLS shape functions and derivatives.
3. In MATLAB, simply run:
   ```matlab
   run('mainMLS1DBar.m')
