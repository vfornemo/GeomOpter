
#![allow(non_upper_case_globals)]

/// Equilibrium C-C bond length r<sub>0</sub> in Angstrom
pub const R0_CC: f64 = 1.5300;
/// Equilibrium C-H bond length r<sub>0</sub> in Angstrom
pub const R0_CH: f64 = 1.1100;
/// Force constant k<sub>b</sub> for C-C bond
pub const Kb_CC: f64 = 300.0;
/// Force constant k<sub>b</sub> for C-H bond
pub const Kb_CH: f64 = 350.0;
/// Equilibrium angle θ<sub>0</sub> for C-C-C in degrees
pub const A0_CCC: f64 = 109.50;
/// Equilibrium angle θ<sub>0</sub> for X-C-X in degrees, X = H, C
pub const A0_XCX: f64 = 109.50;
/// Force constant k<sub>a</sub> for C-C-C angle
pub const Ka_CCC: f64 = 60.00;
/// Force constant k<sub>a</sub> for X-C-X angle, X = H, C
pub const Ka_XCX: f64 = 35.00;
/// Energy term oscillation frequency for dihedral
pub const N_XCCX: f64 = 3.0;
/// Barrier height for dihedral
pub const A_XCCX: f64 = 0.300;
pub const EPSILON_H: f64 = 0.0300;
pub const EPSILON_C: f64 = 0.0700;
pub const SIGMA_H: f64 = 1.20000;
pub const SIGMA_C: f64 = 1.75000;
/// Wolfe condition parameter
pub const WOLFE_C1: f64 = 0.1;
/// Maximum step size for p<sub>q,k</sub> in internal coordinate optimization
pub const MAX_STEP_SIZE: f64 = 0.02;