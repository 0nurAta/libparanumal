import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# user inputs
# -------------------------------------------------
filename = "bns_surfcoeffs_rank0_2.dat"   # change this
xc, yc, zc = 0.0, 0.0, 0.0               # sphere center

# duplicate merge tolerance
tol = 1e-10

# -------------------------------------------------
# read file
# -------------------------------------------------
data = np.loadtxt(filename, comments="#")

# expected columns:
#  0 time
#  1 elem
#  2 face
#  3 node
#  4 x
#  5 y
#  6 z
#  7 nx
#  8 ny
#  9 nz
# 10 rho
# 11 p
# 12 Cp
# 13 tvx
# 14 tvy
# 15 tvz
# 16 tauwx
# 17 tauwy
# 18 tauwz
# 19 tauw
# 20 Cf

time = data[:, 0]
x    = data[:, 4]
y    = data[:, 5]
z    = data[:, 6]
nx   = data[:, 7]
ny   = data[:, 8]
nz   = data[:, 9]
rho  = data[:,10]
p    = data[:,11]
Cp   = data[:,12]
tvx  = data[:,13]
tvy  = data[:,14]
tvz  = data[:,15]
tauwx = data[:,16]
tauwy = data[:,17]
tauwz = data[:,18]
tauw  = data[:,19]
Cf    = data[:,20]

# -------------------------------------------------
# if multiple times exist in one file, pick first
# -------------------------------------------------
t0 = time[0]
mask = np.isclose(time, t0)

x = x[mask]
y = y[mask]
z = z[mask]
nx = nx[mask]
ny = ny[mask]
nz = nz[mask]
rho = rho[mask]
p = p[mask]
Cp = Cp[mask]
tvx = tvx[mask]
tvy = tvy[mask]
tvz = tvz[mask]
tauwx = tauwx[mask]
tauwy = tauwy[mask]
tauwz = tauwz[mask]
tauw = tauw[mask]
Cf = Cf[mask]

# -------------------------------------------------
# geometric coordinates relative to center
# -------------------------------------------------
xr = x - xc
yr = y - yc
zr = z - zc
r  = np.sqrt(xr*xr + yr*yr + zr*zr)

# azimuth: [-pi, pi]
phi = np.arctan2(yr, xr)

# polar angle from +z: [0, pi]
theta = np.arccos(np.clip(zr / r, -1.0, 1.0))

phi_deg = np.degrees(phi)
theta_deg = np.degrees(theta)

# -------------------------------------------------
# remove duplicate surface points
# -------------------------------------------------
x_round = np.round(x / tol) * tol
y_round = np.round(y / tol) * tol
z_round = np.round(z / tol) * tol

xyz = np.column_stack((x_round, y_round, z_round))
xyz_unique, inverse = np.unique(xyz, axis=0, return_inverse=True)
n_unique = xyz_unique.shape[0]

def grouped_mean(arr, inv, ngrp):
    out = np.zeros(ngrp)
    cnt = np.zeros(ngrp)
    for i, g in enumerate(inv):
        out[g] += arr[i]
        cnt[g] += 1.0
    out /= cnt
    return out

x_u = grouped_mean(x, inverse, n_unique)
y_u = grouped_mean(y, inverse, n_unique)
z_u = grouped_mean(z, inverse, n_unique)

nx_u = grouped_mean(nx, inverse, n_unique)
ny_u = grouped_mean(ny, inverse, n_unique)
nz_u = grouped_mean(nz, inverse, n_unique)

rho_u = grouped_mean(rho, inverse, n_unique)
p_u   = grouped_mean(p, inverse, n_unique)
Cp_u  = grouped_mean(Cp, inverse, n_unique)

tvx_u = grouped_mean(tvx, inverse, n_unique)
tvy_u = grouped_mean(tvy, inverse, n_unique)
tvz_u = grouped_mean(tvz, inverse, n_unique)

tauwx_u = grouped_mean(tauwx, inverse, n_unique)
tauwy_u = grouped_mean(tauwy, inverse, n_unique)
tauwz_u = grouped_mean(tauwz, inverse, n_unique)

tauw_u = grouped_mean(tauw, inverse, n_unique)
Cf_u   = grouped_mean(Cf, inverse, n_unique)

xr_u = x_u - xc
yr_u = y_u - yc
zr_u = z_u - zc
r_u  = np.sqrt(xr_u*xr_u + yr_u*yr_u + zr_u*zr_u)

phi_u = np.arctan2(yr_u, xr_u)
theta_u = np.arccos(np.clip(zr_u / r_u, -1.0, 1.0))

phi_deg_u = np.degrees(phi_u)
theta_deg_u = np.degrees(theta_u)

print(f"time = {t0}")
print(f"raw points    = {len(x)}")
print(f"unique points = {len(x_u)}")
print(f"radius min/max = {np.min(r_u):.12e}, {np.max(r_u):.12e}")

# -------------------------------------------------
# plot Cp versus azimuth for equatorial band
# -------------------------------------------------
band_deg = 10.0
equator_mask = np.abs(theta_deg_u - 90.0) < band_deg

phi_eq = phi_deg_u[equator_mask]
Cp_eq  = Cp_u[equator_mask]
Cf_eq  = Cf_u[equator_mask]

idx = np.argsort(phi_eq)
phi_eq = phi_eq[idx]
Cp_eq  = Cp_eq[idx]
Cf_eq  = Cf_eq[idx]

plt.figure(figsize=(8,5))
plt.plot(phi_eq, Cp_eq, "o-")
plt.xlabel(r"$\phi$ (deg), equatorial band")
plt.ylabel(r"$C_p$")
plt.title(f"Cp on equatorial band at time = {t0}")
plt.grid(True)
plt.tight_layout()

plt.figure(figsize=(8,5))
plt.plot(phi_eq, Cf_eq, "o-")
plt.xlabel(r"$\phi$ (deg), equatorial band")
plt.ylabel(r"$C_f$")
plt.title(f"Cf on equatorial band at time = {t0}")
plt.grid(True)
plt.tight_layout()

# -------------------------------------------------
# 3D scatter colored by Cp
# -------------------------------------------------
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(x_u, y_u, z_u, c=Cp_u, s=8)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Surface points colored by Cp")
ax.set_box_aspect((1,1,1))
fig.colorbar(sc, ax=ax, shrink=0.8, label="Cp")
plt.tight_layout()

# -------------------------------------------------
# 3D scatter colored by Cf
# -------------------------------------------------
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(x_u, y_u, z_u, c=Cf_u, s=8)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Surface points colored by Cf")
ax.set_box_aspect((1,1,1))
fig.colorbar(sc, ax=ax, shrink=0.8, label="Cf")
plt.tight_layout()

plt.show()
