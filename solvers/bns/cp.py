import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# user inputs
# -----------------------------
filename = "bns_surfcoeffs_rank0_30.dat"   # change as needed
xc, yc   = 0.0, 0.0                       # cylinder center
rho_inf  = 1.0
U_inf    = 0.1
c        = 1.0

# duplicate-point tolerance
tol = 1e-10

# -----------------------------
# read file
# -----------------------------
data = np.loadtxt(filename, comments="#")

# expected columns:
# 0 time
# 1 elem
# 2 face
# 3 node
# 4 x
# 5 y
# 6 nx
# 7 ny
# 8 tx
# 9 ty
# 10 rho
# 11 p
# 12 Cp
# 13 tauw
# 14 Cf

time = data[:, 0]
x    = data[:, 4]
y    = data[:, 5]
nx   = data[:, 6]
ny   = data[:, 7]
tx   = data[:, 8]
ty   = data[:, 9]
rho  = data[:,10]
p    = data[:,11]
Cp   = data[:,12]
tauw = data[:,13]
Cf   = data[:,14]

# if multiple times are present, pick the first one
t0 = time[0]
mask = np.isclose(time, t0)

x    = x[mask]
y    = y[mask]
nx   = nx[mask]
ny   = ny[mask]
tx   = tx[mask]
ty   = ty[mask]
rho  = rho[mask]
p    = p[mask]
Cp   = Cp[mask]
tauw = tauw[mask]
Cf   = Cf[mask]

# -----------------------------
# recompute Cp
# -----------------------------
p_inf = rho_inf * c**2
dynp  = 0.5 * rho_inf * U_inf**2
Cp_recomputed = (p - p_inf) / dynp

theta = np.arctan2(y - yc, x - xc)
theta_wrap = np.mod(theta, 2.0*np.pi)
theta_deg = np.degrees(theta_wrap)

r = np.sqrt((x - xc)**2 + (y - yc)**2)

# -----------------------------
# remove duplicate geometric points
# -----------------------------
x_round = np.round(x / tol) * tol
y_round = np.round(y / tol) * tol

xy = np.column_stack((x_round, y_round))
xy_unique, inverse = np.unique(xy, axis=0, return_inverse=True)

n_unique = xy_unique.shape[0]

def grouped_mean(arr, inv, ngrp):
    out = np.zeros(ngrp)
    cnt = np.zeros(ngrp)
    for i, g in enumerate(inv):
        out[g] += arr[i]
        cnt[g] += 1.0
    out /= cnt
    return out

x_u    = grouped_mean(x, inverse, n_unique)
y_u    = grouped_mean(y, inverse, n_unique)
nx_u   = grouped_mean(nx, inverse, n_unique)
ny_u   = grouped_mean(ny, inverse, n_unique)
tx_u   = grouped_mean(tx, inverse, n_unique)
ty_u   = grouped_mean(ty, inverse, n_unique)
rho_u  = grouped_mean(rho, inverse, n_unique)
p_u    = grouped_mean(p, inverse, n_unique)
Cp_u   = grouped_mean(Cp, inverse, n_unique)
Cpr_u  = grouped_mean(Cp_recomputed, inverse, n_unique)
tauw_u = grouped_mean(tauw, inverse, n_unique)
Cf_u   = grouped_mean(Cf, inverse, n_unique)
r_u    = grouped_mean(r, inverse, n_unique)

theta_u = np.arctan2(y_u - yc, x_u - xc)
theta_u = np.mod(theta_u, 2.0*np.pi)
theta_deg_u = np.degrees(theta_u)

# sort by angle
idx = np.argsort(theta_u)

x_u    = x_u[idx]
y_u    = y_u[idx]
nx_u   = nx_u[idx]
ny_u   = ny_u[idx]
tx_u   = tx_u[idx]
ty_u   = ty_u[idx]
rho_u  = rho_u[idx]
p_u    = p_u[idx]
Cp_u   = Cp_u[idx]
Cpr_u  = Cpr_u[idx]
tauw_u = tauw_u[idx]
Cf_u   = Cf_u[idx]
r_u    = r_u[idx]
theta_u = theta_u[idx]
theta_deg_u = theta_deg_u[idx]

# -----------------------------
# quick checks
# -----------------------------
cp_err = np.max(np.abs(Cp_u - Cpr_u))
print(f"time = {t0}")
print(f"number of raw points      = {len(x)}")
print(f"number of unique points   = {len(x_u)}")
print(f"radius min/max            = {np.min(r_u):.12e}, {np.max(r_u):.12e}")
print(f"max |Cp - Cp_recomputed|  = {cp_err:.12e}")

# -----------------------------
# optional potential-flow curve
# -----------------------------
theta_ref = np.linspace(0.0, 2.0*np.pi, 400)
cp_potential = 1.0 - 4.0*np.sin(theta_ref)**2

# -----------------------------
# plots
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(theta_deg_u, Cp_u, "o-", label="stored Cp")
plt.plot(theta_deg_u, Cpr_u, "--", label="recomputed Cp")
plt.plot(np.degrees(theta_ref), cp_potential, "k:", label="potential flow")
plt.xlabel(r"$\theta$ (deg)")
plt.ylabel(r"$C_p$")
plt.title(f"Cp distribution at time = {t0}")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.figure(figsize=(8,5))
plt.plot(theta_deg_u, Cf_u, "o-")
plt.xlabel(r"$\theta$ (deg)")
plt.ylabel(r"$C_f$")
plt.title(f"Cf distribution at time = {t0}")
plt.grid(True)
plt.tight_layout()

plt.figure(figsize=(6,6))
plt.plot(x_u, y_u, "o-")
plt.axis("equal")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Unique surface points used for plotting")
plt.grid(True)
plt.tight_layout()

plt.show()
