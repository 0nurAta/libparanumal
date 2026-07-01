#define ADVECTION_SPEED_X 0.00
#define ADVECTION_SPEED_Y 0.0
#define x0 1.0
#define y0 1.0
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
  if(x>0.4 && x<0.65 && y>1.1 && y<1.5){\
  *(q) = 1000.; \
  } \
  else if(x>0.2 && x<0.85 && y>1.5 && y<1.8){ \
  *(q) = 1000.; \
  } \
  else if(x>1.05 && x<1.3 && y>1.1 && y<1.75){ \
  *(q) = 1000.; \
  }\
  else if(x>1.45 && x<1.70 && y>1.1 && y<1.75){ \
  *(q) = 1000.; \
  }  \
  else if(x>1.3 && x<1.45 && y>1.35 && y<1.55){ \
  *(q) = 1000.; \
  }  else *(q)=0.; \
}
