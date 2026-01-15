"""
Author: Jino George. Group of Dr. Chang Liu

2-D version (x–y only) for the periodic-hill / VPF run in Dedalus.
"""

import numpy as np
import dedalus.public as d3
import logging
import matplotlib.pyplot as plt
from mpi4py import MPI
from scipy.special import erf
import time

logger = logging.getLogger(__name__)
 
# ------------------------- Parameters -------------------------
epsilon = 0.5
Lx, Ly = (9.0, 2.0 + 2*epsilon)

Re = 190
mask_option = "opt_erf_peclet"     # {"erf","opt_erf","opt_erf_new","opt_erf_peclet","shifted","opt_tanh"}
distance_fn = "flexible"           # {"const","flexible","optimum"}

dtype = np.float64
timestepper = d3.RK443
max_timestep = 0.005/2
stop_sim_time = 10

# Resolution (2-D)
#nx, ny = 352, 512
nx, ny = 384, 576
# Interface params
c = 1.0
mask_const = Re / c  # 1/eta

# Restart control
restart = 0
checkpoint_path = '/mnt/e/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re100_nek5000_hill_1c_checkpoint/checkpoints/checkpoints_s1.h5'
load_time = 99

# --------------------- Coordinates, bases, domain ---------------------
coords = d3.CartesianCoordinates('x', 'y')
dist   = d3.Distributor(coords, dtype=np.float64)

xbasis = d3.RealFourier(coords['x'], size=nx, bounds=(0, Lx), dealias=3/2) 
ybasis = d3.Chebyshev(  coords['y'], size=ny, bounds=(-Ly/2, Ly/2), dealias=3/2)

print(f"xbasis size = {xbasis.size}, ybasis size = {ybasis.size}")
print("Re = ", Re)

# ---------------------------- Fields ---------------------------------
p   = dist.Field(name='p',   bases=(xbasis, ybasis))
u   = dist.VectorField(coords, name='u', bases=(xbasis, ybasis))

# tau fields live on boundary (no y-basis); in 2-D we keep x-basis only
tau_u1 = dist.VectorField(coords, name='tau_u1', bases=xbasis)
tau_u2 = dist.VectorField(coords, name='tau_u2', bases=xbasis)
tau_p  = dist.Field(name='tau_p')

mask   = dist.Field(name='mask', bases=(xbasis, ybasis))
inv_k  = dist.Field(name='inv_k', bases=(xbasis, ybasis))  # if needed

dPdx = dist.Field(name='dPdx')     # constant field
dPdx.change_scales(1)   # allocate grid

dPdx['g'] = -2/Re
print(dPdx['g'].shape)
# --------------------------- Local grids ------------------------------
# Shapes are (nx_local, ny)
x, y = dist.local_grids(xbasis, ybasis)

# Left edges
global_left_x = xbasis.bounds[0]
global_left_y = ybasis.bounds[0]

nx_global = xbasis.shape[0]
ny_global = ybasis.shape[0]
dy_val = (Ly)/ny_global

if distance_fn == "const":
    mask_threshold = 6*0.5*dy_val  # Adjust as needed, change from 0.5dy to 6*0.5 * dy for the last case
    
elif distance_fn == "flexible":
    # After you define y_local = np.squeeze(y[0, :, 0])  # shape (ny_local,)
    
    dx = Lx/nx
    y_line = y[0, :]                    # shape (ny_local,)
    ny_loc = y_line.shape[0]

    if ny_loc >= 2:       
        dy_row = np.empty_like(y_line)
        dy_row[1:-1] = 0.5 * (y_line[2:] - y_line[:-2])
        dy_row[0]    = (y_line[1]   - y_line[0])
        dy_row[-1]   = (y_line[-1]  - y_line[-2])
        local_max_dy = float(np.max(np.abs(dy_row)))
    else:
        # ny_local is 0 or 1 on this rank → no meaningful dy; contribute 0 to the global max
        local_max_dy = 0.0
    
    comm = MPI.COMM_WORLD
    #mask_threshold = 2.0 * local_max_dy    # use rank-0's local value to match old behavior
    #mask_threshold = MPI.COMM_WORLD.bcast(mask_threshold, root=0)
    global_max_dy = MPI.COMM_WORLD.allreduce(local_max_dy, op=MPI.MAX)
    mask_threshold = 2.0 * global_max_dy
    
    print("mask_threshold_flexible =", mask_threshold)
    
elif distance_fn == "optimum":
         # Compute damping length and optimal smoothing
    eta       = c*(1.0/Re)
    #gamma      = np.sqrt((1.0 / Re) * (1.0 / mask_const_ref))
    gamma      = np.sqrt((eta/ Re))
    # consider gamma = 0.01 for Re = 100 as a standard number ---> delta_star = 3.11346786*0.01
    delta_star = 3.11346786 * gamma    
    mask_threshold = 2.0 * delta_star      # scalar; auto-broadcasts over d_perp
    print("mask_threshold_optimum =", mask_threshold)
    #mask_threshold = 2.0 * (3.11346786*0.01)
    
##mask_threshold = 1/10*(8*0.5*dy)  # This comes out to 0.0078. The minimum grid spacing esp closer to the wall is 0.000272. so this covers 
## approximately 10 or 11 cells as part of the interface. this is much better than before where 100 cells were covered. the challenge lies
## in the fact that this will make stabilizing the code more difficult. 

# Constructing the global y range using the bounds and local values
global_y = np.linspace(global_left_y, global_left_y + Ly, ny_global)
print("nx_global = ", nx_global)
print("ny_global = ", ny_global)
print(global_left_x)
print(global_left_y)

# Helpers for operators
ex, ey = coords.unit_vector_fields(dist)
lift_basis = ybasis.derivative_basis(1)
lift   = lambda A: d3.Lift(A, lift_basis, -1)
grad_u = d3.grad(u) - ey*lift(tau_u1)

x_average = lambda A: d3.Average(A,'x')

dx = lambda A: d3.Differentiate(A, coords['x'])
dy = lambda A: d3.Differentiate(A, coords['y'])

dudx = dx(u@ex); dudy = dy(u@ex)
dvdx = dx(u@ey); dvdy = dy(u@ey)
omega = dvdx - dudy

avg_mass_flow_rate = d3.integ(u@ex)

# ---------------------- Geometry / distance fn -----------------------
y0 = 1.0
A1 = epsilon

A = 4.5
B = 3.5
C = 1.0/6.0      # From the notes
##y_tanh = -y0 + A1 * np.tanh(3.5*np.abs(x-4.5)- 3.5)  # First wave
y_tanh = 3*C * (1.0 + np.tanh(B*(np.abs(x - A) - B)))


# BUlk velocity, Area of fluid and flowrate computation 
Ubulk = 1.0
#Area_fluid = 25.484925
Area_fluid = 25.48

Hbar = (Area_fluid/ Lx)
Q = float(Ubulk * Area_fluid)
beta = 2

print("area-fluid = ", Area_fluid)
print("average height = ",Hbar)
print("flowrate = ",Q)
# --- helpers for signed distance ---
def _wrap(s, L): return s % L
def _pdelta(a, b, L): return ((a - b + 0.5*L) % L) - 0.5*L

def _yb(s, y_base, A, Lx):
    s = _wrap(s, Lx)
    a = 3.5 * (np.abs(s - 0.5*Lx) - 3.5)
    return y_base + A * np.tanh(a)

def _dyb(s, y_base, A, Lx):
    s = _wrap(s, Lx)
    a = 3.5 * (np.abs(s - 0.5*Lx) - 3.5)
    sech2 = 1.0 / np.cosh(a)**2
    sign = np.where(s >= 0.5*Lx, 1.0, -1.0)
    return A * (3.5 * sign) * sech2

def _ddyb(s, y_base, A, Lx):
    s = _wrap(s, Lx)
    a = 3.5 * (np.abs(s - 0.5*Lx) - 3.5)
    sech2 = 1.0 / np.cosh(a)**2
    ap = 3.5 * np.where(s >= 0.5*Lx, 1.0, -1.0)
    return A * (-2.0) * sech2 * np.tanh(a) * (ap**2)

def signed_distance_to_wavy_walls_sampled(x_pt, y_pt, y_base1, A1, Lx, nsample=None):
    s = _wrap(x_pt, Lx)
    for _ in range(8):
        yb  = _yb(s,  y_base1, A1, Lx)
        yp  = _dyb(s, y_base1, A1, Lx)
        ypp = _ddyb(s, y_base1, A1, Lx)
        f   = _pdelta(s, x_pt, Lx) - (y_pt - yb) * yp
        fp  = 1.0 + yp*yp - (y_pt - yb) * ypp
        fp  = np.where(np.abs(fp) < 1e-12, 1e-12, fp)
        s   = _wrap(s - 0.7*(f/fp), Lx)

    yb = _yb(s, y_base1, A1, Lx)
    dx = _pdelta(s, x_pt, Lx)
    dy = y_pt - yb
    d  = np.sqrt(dx*dx + dy*dy)
    phi = y_pt - _yb(_wrap(x_pt, Lx), y_base1, A1, Lx)
    return np.where(phi > 0.0, d, -d)



# ------------------- d_perp & mask (LOCAL, 2-D) ----------------------
## d_perp parallelized code 

y_local = np.squeeze(y[0, :])       # shape (ny_local,)

d_perp = np.zeros((len(y_local), nx))  # shape = (ny_local, nx)

for j in range(len(y_local)):
    for i in range(nx):
        x_pt = x[i]  # Only need x/y from mid-z plane
        y_pt = y_local[j]
        d_perp[j, i] = signed_distance_to_wavy_walls_sampled(x_pt, y_pt, -y0, A1, Lx, nsample=None)

print("Computed d_perp values successfully.")


# Compute binary mask for symmetric wavy walls

# Define binary region identifiers (not normalized values)
solid_region = d_perp < -mask_threshold
fluid_region = d_perp > mask_threshold
interface_region = ~np.logical_or(solid_region, fluid_region)

# Optional visualization (triple mask values: 1 for solid, 0.5 for interface, 0 for fluid)
binary_mask = np.zeros_like(d_perp)
binary_mask[solid_region] = 1
binary_mask[interface_region] = 0.5
binary_mask[fluid_region] = 0

print("Binary mask computed successfully.")

comm = MPI.COMM_WORLD
rank = comm.rank

# --- Gather mask data on rank 0 ---
gathered_mask = comm.gather(binary_mask, root=0)

# These are in shape (nx_local, ny_local), so we transpose them
gathered_x = comm.gather(np.transpose(x[:, :]), root=0)  # shape (ny_local, nx)
gathered_y = comm.gather(np.transpose(y[:, :]), root=0)  # shape (ny_local, nx)

if rank == 0:
    full_mask = np.concatenate(gathered_mask, axis=0)
    full_y = np.concatenate(gathered_y, axis=0)
    full_x = gathered_x[0]
    full_x, full_y = np.meshgrid(full_x.flatten(), full_y.flatten())

    print("full_mask shape:", full_mask.shape)
    print("full_x shape:", full_x.shape)
    print("full_y shape:", full_y.shape)
    shifted_y1 = full_y + 1.5
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(full_x, shifted_y1, full_mask, shading='auto', cmap='gray_r')
    cbar = plt.colorbar(label="Mask Value")
    cbar.set_label("Mask Value", fontsize=22)  
    cbar.ax.tick_params(labelsize=20)  # Increase tick label font size
    plt.xlabel("x",fontsize=24)
    plt.ylabel("y",fontsize=24)
    plt.xticks(fontsize=24)
    plt.yticks([0, 1, 2, 3], fontsize=24)
    #plt.title("Binary Mask for Symmetric Wavy Walls (MPI)")
    plt.tight_layout()
    plt.grid(False)
    plt.savefig("binary_mask_preview.png", dpi=150)
    print("Saved plot to binary_mask_preview.png")  

    try:
        print("Paused: Press Enter to continue...")
    except EOFError:
        pass  # Avoid crash if running without stdin (e.g. mpiexec with no terminal)



def normalized_mask1(d_perp, mask_threshold):
    return 0.5 * (1 - erf(np.sqrt(np.pi) * d_perp / mask_threshold)) * mask_const

def normalized_mask1_shifted(d_perp, mask_threshold, mask_const):
    shift = mask_threshold / np.sqrt(np.pi)
    d_s   = d_perp + shift
    return 0.5 * (1 - erf(np.sqrt(np.pi) * d_s / mask_threshold)) * mask_const

def normalized_mask_optimized(d_perp, Re, mask_const, c):
    eta        = c*(1.0/Re)
    gamma      = np.sqrt(eta/ Re)
    delta_star = 3.11346786 * gamma
    return 0.5 * (1 - erf(np.sqrt(np.pi) * d_perp / delta_star)) * mask_const

def normalized_mask_optimized2(d_perp, Re, mask_const, c):
    j = np.arange(ny_global)
    y_cheb  = np.cos(np.pi * j / (ny_global - 1))
    Y2      = (Ly/2) * y_cheb
    dy_wall = abs(Y2[1] - Y2[0])
    dx_wall = Lx / nx
    h_wall  = min(dx_wall, dy_wall)
    k_cells = 6.0
    delta_floor = k_cells * h_wall

    delta_opt = 3.11346786 * np.sqrt(c) / Re

    s_min, Re_t, p = 0.6, 40.0, 2.0
    s = s_min + (1.0 - s_min) / (1.0 + (Re_t / Re)**p)
    if Re > 80:
        s = 1 + (1.0 - s_min) / (1.0 + (Re_t / Re)**p)

    delta_star = max(s * delta_opt, delta_floor)
    return 0.5 * (1 - erf(np.sqrt(np.pi) * d_perp / delta_star)) * mask_const

def normalized_mask_optimized3(d_perp, Re, mask_const, c,
                               Hav=1.0, Pe_max=2.0, kappa=1.0):
    j = np.arange(ny_global)
    y_cheb  = np.cos(np.pi * j / (ny_global - 1))
    Y2      = (Ly/2) * y_cheb
    dy_wall = abs(Y2[1] - Y2[0])
    dx_wall = Lx / nx
    h_wall  = min(dx_wall, dy_wall)
    k_cells = 50.0
    delta_floor = k_cells * h_wall

    delta_opt = 3.11346786 * (np.sqrt(c) / Re)

    s_min, Re_t, p = 0.6, 80.0, 2.0
    s_lo  = s_min + (1.0 - s_min) / (1.0 + (Re_t / Re)**p)
    s     = 1.0 if Re > 45 else s_lo

    delta_burns = s * delta_opt
    # optional Pe cap was disabled in your version
    delta_star  = max(delta_burns, delta_floor)

    U_t_est  = kappa / Hav
    Pe_delta = U_t_est * delta_star * Re
    print("delta_floor=", delta_floor, " delta_opt=", delta_opt,
          " delta_burns=", delta_burns, " delta_star=", delta_star,
          " Pe_delta=", Pe_delta,
          " s=", s)
    
    return 0.5 * (1 - erf(np.sqrt(np.pi) * d_perp / delta_star)) * mask_const

def normalized_mask_optimized_tanh(d_perp, Re, mask_const, c):
    eta = 1.0 / (Re / c)
    gamma = np.sqrt(eta / Re)
    delta_star = 2.64822828 * gamma
    return 0.5 * (1 - np.tanh(d_perp / delta_star)) * mask_const

# Choose mask
if mask_option == "erf":
    final_mask_2d_local = normalized_mask1(d_perp, mask_threshold)
elif mask_option == "opt_erf":
    final_mask_2d_local = normalized_mask_optimized(d_perp, Re, mask_const, c)
elif mask_option == "opt_erf_new":
    final_mask_2d_local = normalized_mask_optimized2(d_perp, Re, mask_const, c)
elif mask_option == "opt_erf_peclet":
    final_mask_2d_local = normalized_mask_optimized3(d_perp, Re, mask_const, c, Hav=1.0, Pe_max=2.0, kappa=1.0)
elif mask_option == "shifted":
    final_mask_2d_local = normalized_mask1_shifted(d_perp, mask_threshold, mask_const)
elif mask_option == "opt_tanh":
    final_mask_2d_local = normalized_mask_optimized_tanh(d_perp, Re, mask_const, c)
else:
    raise ValueError(f"Unknown mask_option='{mask_option}'")


print(final_mask_2d_local.shape)

comm = MPI.COMM_WORLD
rank = comm.rank

gathered_mask_norm = comm.gather(final_mask_2d_local, root=0)
gathered_x = comm.gather(np.transpose(x[:, :]), root=0)  # shape (ny_local, nx)
gathered_y = comm.gather(np.transpose(y[:, :]), root=0)  # shape (ny_local, nx)

if rank == 0:
    print("This is rank 0 â€” safe to do plotting or file I/O here.")
    print("gathered_mask_norm length:", len(gathered_mask_norm))
    print("gathered_x[0] shape:", gathered_x[0].shape)
    print("gathered_y[0] shape:", gathered_y[0].shape)
    

if rank == 0:
    full_mask = np.concatenate(gathered_mask_norm, axis=0)  # shape: (ny_total, nx)
    full_y = np.concatenate(gathered_y, axis=0)              # shape: (ny_total, nx)
    full_x = gathered_x[0]                                   # all ranks have same x

    # Create meshgrid for plotting
    full_x, full_y = np.meshgrid(full_x.flatten(), full_y.flatten())

    print("full_mask shape:", full_mask.shape)
    print("full_x shape:", full_x.shape)
    print("full_y shape:", full_y.shape)
    print("full mask max:", np.max(full_mask))
    print("full mask min:", np.min(full_mask))

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import matplotlib.ticker as ticker


    #cmap = ListedColormap(['#ADD8E6', 'black'])
    cmap = ListedColormap(['#ADD8E6', '#FF0000', 'black'])  # light blue, black, red

    plt.figure(figsize=(8, 6))
    shifted_y = full_y + 1.5

    plt.pcolormesh(full_x, shifted_y, full_mask, shading='auto', cmap='gray_r')
    #plt.pcolormesh(full_x, shifted_y, full_mask, shading='auto', cmap=cmap)
    #plt.pcolormesh(full_x, shifted_y, full_mask, shading='auto', cmap='viridis')

    #plt.colorbar(label="Mask Value")
    
    cbar = plt.colorbar(label="Mask Value")
    cbar.set_label("Mask Value", fontsize=22)  
    #cbar.set_ticks([0, 100])           # only two ticks
    #cbar.set_ticklabels(["0", "100"])  # corresponding labels
    cbar.ax.yaxis.set_major_formatter(
    ticker.FuncFormatter(lambda val, pos: f"{val/100:.1f}")
    )

    cbar.ax.tick_params(labelsize=20)  # Increase tick label font size
    plt.xlabel("x",fontsize=24)
    plt.ylabel("y",fontsize=24)
    #plt.xticks(fontsize=24)
    plt.xticks([0, 3, 6, 9], fontsize=24)
    plt.yticks([0, 1, 2, 3], fontsize=24)
    #plt.title("Normalized Mask for Solid/Fluid Interface (MPI)")
    plt.tight_layout()
    plt.grid(False)
    plt.savefig("normalized_mask_preview.png", dpi=150)
    print("âœ… Saved plot to normalized_mask_preview.png")

    time.sleep(2)  # Optional: short pause


# Assign to Dedalus field:
# final_mask_2d_local is (ny, nx_local) -> Dedalus expects (nx_local, ny)
mask['g'] = final_mask_2d_local.T

# Sanity checks across ranks
comm  = MPI.COMM_WORLD
rank  = comm.rank

if mask['g'].size > 0:
    local_min = np.min(mask['g']); local_max = np.max(mask['g'])
else:
    local_min, local_max = np.inf, -np.inf

global_min = comm.allreduce(local_min, op=MPI.MIN)
global_max = comm.allreduce(local_max, op=MPI.MAX)
if rank == 0:
    print("✅ Mask field assigned. Global range:", global_min, global_max)

# --------------------------- Problem ---------------------------------
problem = d3.IVP([p, u, dPdx, tau_p, tau_u1, tau_u2], namespace=locals())
problem.add_equation("trace(grad_u) + tau_p = 0")
problem.add_equation("dt(u) - 1/Re*div(grad_u) + grad(p) + lift(tau_u2) + dPdx*ex = -dot(u,grad(u)) - mask*u")
problem.add_equation("u(y=-Ly/2) = 0")
problem.add_equation("u(y=+Ly/2) = 0")
problem.add_equation("integ(p) = 0")
problem.add_equation("integ(u@ex) = Q")   # keep your flow-rate constraint

print("Number of equations:", len(problem.equations))
print("Number of variables:", len(problem.variables))

# -------------------------- Solver & IO -------------------------------
dt = 1e-4
solver = problem.build_solver(timestepper)
solver.stop_sim_time = 1000

# Restart / ICs
if restart:
    wrote, initial_dt = solver.load_state(checkpoint_path, allow_missing=True)
    file_handler_mode = 'append'
else:
    amplitude = 1e-6
    #u['g'][0] = (1 - y**2) + np.random.randn(*u['g'][0].shape) * amplitude * np.sin(np.pi*(y+Ly/2)/Ly)
    u['g'][0] = 1 - (y/(Ly/2))**2 + np.random.randn(*u['g'][0].shape) * amplitude * np.sin(np.pi*(y+Ly/2)/Ly)
    u['g'][1] = 0.0
    dPdx['g'] = -2/Re
    file_handler_mode = 'overwrite'

# Checkpoints (full state)
checkpoints = solver.evaluator.add_file_handler('checkpoints', iter=2000, max_writes=999)
for fld in solver.state:
    checkpoints.add_task(fld, name=fld.name, layout='g')

# Snapshots (2-D tasks only)
snaps = solver.evaluator.add_file_handler('snapshots_channel', sim_dt=1, max_writes=400)
snaps.add_task(u,     name='velocity')
snaps.add_task(dPdx,  name='dPdx')
snaps.add_task(dudx,  name='dUdx')
snaps.add_task(dudy,  name='dUdy')
snaps.add_task(dvdx,  name='dVdx')
snaps.add_task(dvdy,  name='dVdy')
snaps.add_task(mask,  name='mask')
snaps.add_task(avg_mass_flow_rate, name='avg_mass_flow_rate')

snaps_stress = solver.evaluator.add_file_handler('snapshots_channel_stress', sim_dt=1, max_writes=400)
snaps_stress.add_task(x_average(u), name='ubar')
snaps_stress.add_task(omega,        name='omega')  # scalar vorticity

# CFL
CFL = d3.CFL(solver, initial_dt=dt, cadence=10, safety=0.3, threshold=0.02,
             max_change=1.1, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties (2-D)
flow = d3.GlobalFlowProperty(solver, cadence=20)
flow.add_property(((u - x_average(u))@(u - x_average(u)))/2, name='TKE')
flow.add_property(x_average((u@u)/2), name='KE')
flow.add_property(d3.integ(u@u/2), name='KE_int')
flow.add_property(
    d3.integ(((u - x_average(u))@(u - x_average(u)))/2),
    name='TKE_int')

flow.add_property(d3.integ(u@ex), name='Q_now')  # volumetric flow rate. this should be constant
# --------------------------- Main loop --------------------------------
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)
        if (solver.iteration-1) % 10 == 0:
            max_TKE = flow.max('TKE')
            KE      = flow.max('KE')
            KE_int  = flow.max('KE_int')
            TKE_int = flow.max('TKE_int')
            Q_now   = flow.max('Q_now')
            #dPdx_val = float(np.mean(dPdx['g']))
            logger.info('Iteration=%i, Time=%e, dt=%e, max(TKE)=%f, KE = %f, KE_int=%f, TKE_int=%f, Q_now = %f',
                        solver.iteration, solver.sim_time, timestep, max_TKE, KE, KE_int, TKE_int, Q_now)
except Exception as e:
    import traceback
    logger.error(f'Exception raised: {e}')
    logger.error(traceback.format_exc())
    raise
finally:
    solver.log_stats()
