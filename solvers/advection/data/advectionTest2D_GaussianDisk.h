#define ADVECTION_SPEED_X 0.05
#define ADVECTION_SPEED_Y 0.0
#define x0 1.0
#define y0 1.0
#define R0 0.40
// Flux function
#define advectionFlux2D(t, x, y, q, cx, cy) \
{                                       \
  *(cx) = ADVECTION_SPEED_X*q;          \
  *(cy) = ADVECTION_SPEED_Y*q;          \
}

// max wavespeed (should be max eigen of Jacobian of flux function)
#define advectionMaxWaveSpeed2D(t, x, y, q, u, v) \
{                                                 \
  *(u) = ADVECTION_SPEED_X;                       \
  *(v) = ADVECTION_SPEED_Y;                       \
}

// Boundary conditions
/* wall 1, outflow 2 */
#define advectionDirichletConditions2D(bc, t, x, y, nx, ny, qM, qB) \
{                                       \
  if(bc==1){                            \
    *(qB) = 0.0;                        \
  } else if(bc==2){                     \
    *(qB) = qM;                         \
  }                                     \
}

// Initial conditions
#define advectionInitialConditions2D(t, x, y, q) \
{                                       \
  *(q) = 2*exp(-80*(sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - R0) \
              *(sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - R0)); \
}
