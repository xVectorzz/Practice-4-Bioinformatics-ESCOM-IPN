import cv2
import numpy as np
import math
import periodictable as pt

# --- Display and Setup Variables ---
WINDOW_NAME = "Simulacion de Fuerzas 2D"
DISPLAY_WAIT_MS = 20 

# Colores (BGR)
BLUE_H = (0, 0, 255) 
RED_O = (255, 0, 0) 
NA_COLOR = (0, 100, 255)

maxX = 500
maxY = 500
x0 = int(maxX / 2)
y0 = int(maxY / 2)

# --- Critical Stability Constants ---
F_MAX = 70.0 
DAMPING = 0.985 
COULOMB_K = 4.0 
R_MIN_COULOMB = 0.5
K_BOND = 50.0 
K_ANGLE = 20.0
SIGMA_LJ = 1.0 
EPSILON_LJ = 1.0 
G_ACCELERATION = 0.05 

# --- Physical and Geometric Constants ---
mass_O = pt.O.mass
mass_H = pt.H.mass
mass_NA = pt.Na.mass
CHARGE_O = -0.84
CHARGE_H = 0.42
CHARGE_NA = 1.0
BOND_LENGTH = 1.0 
BOND_ANGLE = np.radians(104.5)
SCALE_FACTOR = 25 
delta_t = 0.03 

# --- Classes ---
class Atom:
    # Represents a single particle with mass and charge.
    def __init__(self, name, mass, charge, r, color, v=np.array([0.0, 0.0])):
        self.name = name
        self.mass = mass
        self.charge = charge 
        self.color = color
        self.r = r 
        self.v = v 
        self.a = np.array([0.0, 0.0])
        self.f = np.array([0.0, 0.0]) 
    def apply_force(self, force):
        self.f += force
    def reset_force(self):
        self.f = np.array([0.0, 0.0])

class Molecule:
    # Represents a collection of atoms with defined bonds and angles.
    def __init__(self, atoms, bonds=None, angles=None):
        self.atoms = atoms
        self.bonds = bonds if bonds is not None else [] 
        self.angles = angles if angles is not None else [] 
    def get_atom(self, index):
        return self.atoms[index]

# --- Utility and Drawing Functions ---

def draw_circle(x, y, color):
    # Draws a circle on the image, adjusting coordinates.
    cv2.circle(img, (int(x0 + x), int(y0 - y)), 10, color, -1)

def create_h2o(center_r):
    # Creates an H2O molecule structure.
    angle = BOND_ANGLE / 2
    r_h1 = center_r + BOND_LENGTH * np.array([np.cos(angle), np.sin(angle)])
    r_h2 = center_r + BOND_LENGTH * np.array([np.cos(angle), -np.sin(angle)])
    
    atom_o = Atom('O', mass_O, CHARGE_O, center_r, RED_O)
    atom_h1 = Atom('H', mass_H, CHARGE_H, r_h1, BLUE_H)
    atom_h2 = Atom('H', mass_H, CHARGE_H, r_h2, BLUE_H)
    
    bonds = [(0, 1, BOND_LENGTH), (0, 2, BOND_LENGTH)]
    angles = [(1, 0, 2, BOND_ANGLE)]
    return Molecule([atom_o, atom_h1, atom_h2], bonds, angles)

def is_intramolecular(atom_i, atom_j, mol_list):
    # Checks if two atoms are in the same molecule.
    mol_i = None
    mol_j = None
    
    for mol in mol_list:
        if atom_i in mol.atoms: mol_i = mol
        if atom_j in mol.atoms: mol_j = mol
    
    if mol_i is not mol_j or mol_i is None:
        return False, None 
        
    return True, mol_i

# --- Force Functions (Core Physics) ---

def hooke_force(r_vec, k, r, r0):
    # Calculates Hooke's spring force (Bonding Force).
    if r == 0: return np.zeros(2)
    direction = r_vec / r
    F = -k * (r - r0) * direction
    return F

def angle_force(p0, pH1, pH2, kTheta=K_ANGLE, theta0=BOND_ANGLE):
    # Calculates the angular force (H-O-H). Returns F0, F1, F2.
    v1 = pH1.r - p0.r
    v2 = pH2.r - p0.r
    r1 = np.linalg.norm(v1)
    r2 = np.linalg.norm(v2)
    if r1==0 or r2==0: return np.zeros(2), np.zeros(2), np.zeros(2) 
    
    cos_theta = np.dot(v1, v2) / (r1 * r2)
    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    dTheta = theta - theta0
    fMag = -kTheta * dTheta 
    
    n1 = np.array([-v1[1], v1[0]]) / (r1 + 1e-9)
    n2 = np.array([ v2[1],-v2[0]]) / (r2 + 1e-9)
    
    F1 = fMag * n1
    F2 = fMag * n2
    F0 = -(F1 + F2) 
    
    return F0, F1, F2 

def lennard_jones_force(p1, p2, epsilon=EPSILON_LJ, sigma=SIGMA_LJ):
    # Calculates Lennard-Jones force (Non-bonded interaction).
    r_vec = p1.r - p2.r
    r = np.linalg.norm(r_vec)
    if r == 0:
        return np.zeros(2), np.zeros(2)
    direction = r_vec / r
    r7 = r**7
    r13 = r**13
    if r7 == 0 or r13 == 0:
        return np.zeros(2), np.zeros(2)
    fMag = 24 * epsilon * ((2 * (sigma**12) / r13) - ((sigma**6) / r7))
    fMag = np.clip(fMag, -F_MAX, F_MAX)
    F = fMag * direction
    return F, -F


def coulomb_force(p1, p2, k=COULOMB_K, minDist=R_MIN_COULOMB):
    # Calculates the Coulomb force (Electrostatic interaction).
    r_vec = p1.r - p2.r 
    r = np.linalg.norm(r_vec)
    if r < minDist:
        r = minDist
    direction = r_vec / r
    F_mag = k * (p1.charge * p2.charge) / (r**2)
    F = F_mag * direction
    F = np.clip(F, -F_MAX, F_MAX)
    return F, -F


# --- Constraints and Step (Simulation) ---

def link_constraint(mol_list, dt):
    # Applies geometric link constraint (maintains O-H bond length).
    for mol in mol_list:
        if mol.atoms[0].name != 'O': continue 

        for i in range(1, 3): 
            p1 = mol.get_atom(0) # Oxygen
            p2 = mol.get_atom(i) # Hydrogen
            
            delta = p2.r - p1.r
            dist = np.linalg.norm(delta)
            if dist == 0: continue
                
            diff = (dist - BOND_LENGTH) / dist
            
            p2.r -= delta * diff * 0.5 
            p1.r += delta * diff * 0.01 


def step(dt=delta_t, damping=DAMPING):
    # Main function to advance the simulation one time step (integration).
    
    for p in all_atoms:
        p.reset_force()
        p.a = np.zeros(2)
        
        # Apply Gravity (currently near zero)
        F_gravity = np.array([0.0, -G_ACCELERATION * p.mass])
        p.apply_force(F_gravity)
    
    water_molecules = [water1, water2]
    
    # STEP 1: Intra-molecular Forces (Bonds and Angles)
    for mol in water_molecules:
        p0, p1, p2 = mol.atoms[0], mol.atoms[1], mol.atoms[2]
        
        # a. Angular Force (H-O-H)
        f0, f1, f2 = angle_force(p0, p1, p2)
        p0.apply_force(f0); p1.apply_force(f1); p2.apply_force(f2)
        
        # b. Bonding Force (Hooke: O-H1 and O-H2)
        r_vec_oh1 = p1.r - p0.r
        r_oh1 = np.linalg.norm(r_vec_oh1)
        f_oh1 = hooke_force(r_vec_oh1, k=K_BOND, r=r_oh1, r0=BOND_LENGTH)
        p0.apply_force(-f_oh1); p1.apply_force(f_oh1)
        
        r_vec_oh2 = p2.r - p0.r
        r_oh2 = np.linalg.norm(r_vec_oh2)
        f_oh2 = hooke_force(r_vec_oh2, k=K_BOND, r=r_oh2, r0=BOND_LENGTH)
        p0.apply_force(-f_oh2); p2.apply_force(f_oh2)


    # STEP 2: Non-bonded Forces (Inter-atomic)
    num_atoms = len(all_atoms)
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            atom_i = all_atoms[i]
            atom_j = all_atoms[j]

            # Check if atoms are directly bonded (skip non-bonded forces if bonded)
            is_intra, mol_container = is_intramolecular(atom_i, atom_j, water_molecules)
            bonded_pair = False

            if is_intra and mol_container is not None:
                for b in mol_container.bonds:
                    if (mol_container.atoms.index(atom_i), mol_container.atoms.index(atom_j)) in [(b[0], b[1]), (b[1], b[0])]:
                        bonded_pair = True
                        break

            if not bonded_pair:
                # Coulomb Force
                f_coul_i, f_coul_j = coulomb_force(atom_i, atom_j)
                atom_i.apply_force(f_coul_i)
                atom_j.apply_force(f_coul_j)

                # Lennard-Jones Force
                f_lj_i, f_lj_j = lennard_jones_force(atom_i, atom_j)
                atom_i.apply_force(f_lj_i)
                atom_j.apply_force(f_lj_j)

                
    # STEP 3: Kinematic Integration
    for p in all_atoms:
        if p.mass > 0:
            p.a = p.f / p.mass

        # Integration (Euler/Verlet) and Damping
        p.v += p.a * dt
        p.v *= damping 
        p.r += p.v * dt 

        # Wall Logic (Boundary reflection)
        max_limit_x = maxX / SCALE_FACTOR / 2
        max_limit_y = maxY / SCALE_FACTOR / 2
        
        # X-reflection
        if p.r[0] > max_limit_x: 
            p.r[0] = max_limit_x
            p.v[0] *= -1.0 
        elif p.r[0] < -max_limit_x: 
            p.r[0] = -max_limit_x
            p.v[0] *= -1.0

        # Y-reflection
        if p.r[1] > max_limit_y: 
            p.r[1] = max_limit_y
            p.v[1] *= -1.0
        elif p.r[1] < -max_limit_y: 
            p.r[1] = -max_limit_y
            p.v[1] *= -1.0

    # STEP 4: Constraints (Geometric stability)
    link_constraint(water_molecules, dt)

def draw():
    # Draws all atoms and their bonds on the OpenCV canvas.
    for atom in all_atoms:
        draw_circle(atom.r[0]*SCALE_FACTOR, atom.r[1]*SCALE_FACTOR, atom.color) 

    def draw_molecule_lines(mol):
        for idx1, idx2, _ in mol.bonds:
            r1 = mol.get_atom(idx1).r
            r2 = mol.get_atom(idx2).r
            cv2.line(img, (int(x0 + r1[0] * SCALE_FACTOR), int(y0 - r1[1] * SCALE_FACTOR)), 
                     (int(x0 + r2[0] * SCALE_FACTOR), int(y0 - r2[1] * SCALE_FACTOR)), 
                     (128, 128, 128), 1)

    draw_molecule_lines(water1)
    draw_molecule_lines(water2)
    
# --- Modified Initialization ---

DIST_FROM_CENTER = 4.0 
V_RADIAL = 1.5 
V_TANGENCIAL = 2.0

# Sodium Ion at the center
r_na = np.array([0.0, 0.0])
sodium_ion = Atom('Na', mass_NA, CHARGE_NA, r_na, NA_COLOR, v=np.array([0.0, 0.0]))

# 1. Water Molecule 1 (Left)
r_o1_center = np.array([-DIST_FROM_CENTER, 0.0])
V_INIT_WATER1 = np.array([V_RADIAL, V_TANGENCIAL]) 
water1 = create_h2o(r_o1_center)

for atom in water1.atoms:
    atom.v = V_INIT_WATER1 

# 2. Water Molecule 2 (Right)
r_o2_center = np.array([DIST_FROM_CENTER, 0.0])
V_INIT_WATER2 = np.array([-V_RADIAL, -V_TANGENCIAL]) 
water2 = create_h2o(r_o2_center)

for atom in water2.atoms:
    atom.v = V_INIT_WATER2 
    # Rotate 180° for molecule 2
    r_atom_rel = atom.r - r_o2_center
    r_rotated_rel = np.array([-r_atom_rel[0], -r_atom_rel[1]])
    atom.r = r_o2_center + r_rotated_rel
    
all_atoms = [sodium_ion] + water1.atoms + water2.atoms 

MaxIterations = 1000 

img = np.zeros((maxX, maxY, 3), dtype="uint8")
cv2.namedWindow(WINDOW_NAME)


# --- Main Simulation and Drawing Loop ---
apply_lj = True 

for count in range(MaxIterations):
    img[:] = (0, 0, 0)

    step()

    draw()

    cv2.imshow(WINDOW_NAME, img)

    key = cv2.waitKey(DISPLAY_WAIT_MS)
    if key == ord('q'):
        break

cv2.destroyAllWindows()