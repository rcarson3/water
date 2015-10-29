#include "stepper.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include <stdio.h>
#include <omp.h>

//ldoc on
/**
 * ## Implementation
 *
 * ### Structure allocation
 */

central2d_t* central2d_init(float w, float h, int nx, int ny,
                            int imin, int imin_, int jmin, int jmin_,
                            int imax, int imax_, int jmax, int jmax_,
                            int nfield, flux_t flux, speed_t speed,
                            float cfl)
{
    // We extend to a four cell buffer to avoid BC comm on odd time steps
    int ng = 4;

    central2d_t* sim = (central2d_t*) malloc(sizeof(central2d_t));
    sim->nx = nx;
    sim->ny = ny;
  
    sim->imin   = imin;
    sim->imin_  = imin_;
    sim->imino  = imin - ng;
    sim->imino_ = imin_ - ng;
  
    sim->jmin   = jmin;
    sim->jmin_  = jmin_;
    sim->jmino  = jmin  - ng;
    sim->jmino_ = jmin_ - ng;
  
    sim->imax   = imax;
    sim->imax_  = imax_;
    sim->imaxo  = imax  + ng;
    sim->imaxo_ = imax_ + ng;
  
    sim->jmax   = jmax;
    sim->jmax_  = jmax_;
    sim->jmaxo  = jmax + ng;
    sim->jmaxo_ = jmax_+ ng;

    sim->ng = ng;
    sim->nfield = nfield;
    sim->dx = w/nx;
    sim->dy = h/ny;
    sim->flux = flux;
    sim->speed = speed;
    sim->cfl = cfl;

    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    int nc = nx_all * ny_all;
    int N  = nfield * nc;
    sim->u  = (float*) malloc((4*N + 6*nx_all)* sizeof(float));
    sim->v  = sim->u +   N;
    sim->f  = sim->u + 2*N;
    sim->g  = sim->u + 3*N;
    sim->scratch = sim->u + 4*N;

    return sim;
}

void central2d_free(central2d_t* sim)
{
    free(sim->u);
    free(sim);
}


int central2d_offset(int nx, int ny, int ng, int k, int ix, int iy)
{
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    return (k*ny_all+(ng+iy))*nx_all+(ng+ix);
}


/**
 * ### Boundary conditions
 *
 * In finite volume methods, boundary conditions are typically applied by
 * setting appropriate values in ghost cells.  For our framework, we will
 * apply periodic boundary conditions; that is, waves that exit one side
 * of the domain will enter from the other side.
 *
 * We apply the conditions by assuming that the cells with coordinates
 * `nghost <= ix <= nx+nghost` and `nghost <= iy <= ny+nghost` are
 * "canonical", and setting the values for all other cells `(ix,iy)`
 * to the corresponding canonical values `(ix+p*nx,iy+q*ny)` for some
 * integers `p` and `q`.
 */

static inline
void copy_subgrid(float* restrict dst,
                  const float* restrict src,
                  int nx, int ny, int stride)
{
  for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix)
      dst[iy*stride+ix] = src[iy*stride+ix];
}

void central2d_pushBC(float* restrict u_, float* restrict u,
                      int nx_, int ny_, int ng,
                      int nx, int ny,
                      int imin_, int jmin_)
{
    printf(" bad stuff 0");
  int nx_all = nx+2*ng;
  int ny_all = ny+2*ng;
  //copy left boundary
  int i = 0;
  for(int k=0;k<=4;k++)
    for(int j=0;j<=ny_; j++){
      u[(k*ny_all+(ng+jmin_+j))*nx_all+(ng+imin_+i)] = u_[(k*ny_all+(ng+j))*nx_all+(ng+i)];
    }

  //copy right boundary
  i = nx_;
  for(int k=0;k<=4;k++)
    for(int j=0;j<=ny_; j++){
      u[(k*ny_all+(ng+jmin_+j))*nx_all+(ng+imin_+i)] = u_[(k*ny_all+(ng+j))*nx_all+(ng+i)];
    }

  //copy bottom boundary
  int j = 0;
    for(int k=0;k<=4;k++)
      for(int i=0;i=nx_; i++){
        u[(k*ny_all+(ng+jmin_+j))*nx_all+(ng+imin_+i)] = u_[(k*ny_all+(ng+j))*nx_all+(ng+i)];
    }

  //copy top boundary
  j = ny_;
  for(int k=0;k<=4;k++)
    for(int i=0;i<=nx_; i++){
      u[(k*ny_all+(ng+jmin_+j))*nx_all+(ng+imin_+i)] = u_[(k*ny_all+(ng+j))*nx_all+(ng+i)];
    }
}

static
void central2d_BCset(float* restrict u_, float* restrict u,
                     int nx_, int ny_,int ng,
                     int nx, int ny,
                     int imin_, int jmin_, int imax_, int jmax_,
                     int imin, int jmin, int imax, int jmax)
{
  //Side BCs
  if(imin_ == imin) { // On left wall, periodic
    int i = 0;
      for(int j=0;j<=ny_;j++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global right side
          u_[central2d_offset(nx_,ny_,ng,0,i,j)-gcell] = u[central2d_offset(nx,ny,ng,0,imax,jmin_+j)-(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,1,i,j)-gcell] = u[central2d_offset(nx,ny,ng,1,imax,jmin_+j)-(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,2,i,j)-gcell] = u[central2d_offset(nx,ny,ng,2,imax,jmin_+j)-(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,3,i,j)-gcell] = u[central2d_offset(nx,ny,ng,3,imax,jmin_+j)-(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,4,i,j)-gcell] = u[central2d_offset(nx,ny,ng,4,imax,jmin_+j)-(gcell-1)];
        }
  }
    else if(imax_ == imax) { // On right wall, periodic
      int i = nx_;
      for(int j=0;j<=ny_;j++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global left side
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+gcell] = u[central2d_offset(nx,ny,ng,0,imin,jmin_+j)+(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,1,i,j)+gcell] = u[central2d_offset(nx,ny,ng,1,imin,jmin_+j)+(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,2,i,j)+gcell] = u[central2d_offset(nx,ny,ng,2,imin,jmin_+j)+(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,3,i,j)+gcell] = u[central2d_offset(nx,ny,ng,3,imin,jmin_+j)+(gcell-1)];
          u_[central2d_offset(nx_,ny_,ng,4,i,j)+gcell] = u[central2d_offset(nx,ny,ng,4,imin,jmin_+j)+(gcell-1)];
        }
    }
    else {//Internal domain
      // Left subdomain side
      int i = 0;
      for(int j=0;j<=ny_;j++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global neighbor
          u_[central2d_offset(nx_,ny_,ng,0,i,j)-gcell] = u[central2d_offset(nx,ny,ng,0,imin_,jmin_+j)-(gcell)];
          u_[central2d_offset(nx_,ny_,ng,1,i,j)-gcell] = u[central2d_offset(nx,ny,ng,1,imin_,jmin_+j)-(gcell)];
          u_[central2d_offset(nx_,ny_,ng,2,i,j)-gcell] = u[central2d_offset(nx,ny,ng,2,imin_,jmin_+j)-(gcell)];
          u_[central2d_offset(nx_,ny_,ng,3,i,j)-gcell] = u[central2d_offset(nx,ny,ng,3,imin_,jmin_+j)-(gcell)];
          u_[central2d_offset(nx_,ny_,ng,4,i,j)-gcell] = u[central2d_offset(nx,ny,ng,4,imin_,jmin_+j)-(gcell)];
        }
    
      // Right subdomain side
      i = nx_;
      for(int j=0;j<=ny_;j++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global neighbor
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+gcell] = u[central2d_offset(nx,ny,ng,0,imax_,jmin_+j)+(gcell)];
          u_[central2d_offset(nx_,ny_,ng,1,i,j)+gcell] = u[central2d_offset(nx,ny,ng,1,imax_,jmin_+j)+(gcell)];
          u_[central2d_offset(nx_,ny_,ng,2,i,j)+gcell] = u[central2d_offset(nx,ny,ng,2,imax_,jmin_+j)+(gcell)];
          u_[central2d_offset(nx_,ny_,ng,3,i,j)+gcell] = u[central2d_offset(nx,ny,ng,3,imax_,jmin_+j)+(gcell)];
          u_[central2d_offset(nx_,ny_,ng,4,i,j)+gcell] = u[central2d_offset(nx,ny,ng,4,imax_,jmin_+j)+(gcell)];
        }
    }
  
    //Top/Bottom BCs
  if(jmin_ == jmin) { // On bottom wall, periodic
    int j = 0;
      for(int i=0;i<=nx_;i++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global top
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax)-(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax)-(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax)-(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax)-(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax)-(gcell-1)*(nx+2*ng)];
        }
  }
    else if(jmax_ == jmax) { // On top wall, periodic
      int j = ny_;
      for(int i=0;i<=nx_;i++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global bottom
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin)+(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin)+(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin)+(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin)+(gcell-1)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin)+(gcell-1)*(nx+2*ng)];
        }
    }
    else {//Internal domain
      // Bottom subdomain side
      int j = 0;
      for(int i=0;i<=nx_;i++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global neighbor
          u_[central2d_offset(nx_,ny_,ng,0,i,j)-(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin_)-(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)-(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin_)-(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)-(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin_)-(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)-(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin_)-(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)-(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmin_)-(gcell)*(nx+2*ng)];
        }
      // Top subdomain side
      j = ny_;
      for(int i=0;i<=nx_;i++)
        for(int gcell=0;gcell<=ng; ++gcell){
          //copy from global neighbor
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax_)+(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax_)+(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax_)+(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax_)+(gcell)*(nx+2*ng)];
          u_[central2d_offset(nx_,ny_,ng,0,i,j)+(gcell*(nx_+2*ng))] = u[central2d_offset(nx,ny,ng,0,imin_+i,jmax_)+(gcell)*(nx+2*ng)];
        }
    
    }
    
}


/**
 * ### Derivatives with limiters
 *
 * In order to advance the time step, we also need to estimate
 * derivatives of the fluxes and the solution values at each cell.
 * In order to maintain stability, we apply a limiter here.
 *
 * The minmod limiter *looks* like it should be expensive to computer,
 * since superficially it seems to require a number of branches.
 * We do something a little tricky, getting rid of the condition
 * on the sign of the arguments using the `copysign` instruction.
 * If the compiler does the "right" thing with `max` and `min`
 * for floating point arguments (translating them to branch-free
 * intrinsic operations), this implementation should be relatively fast.
 */


// Branch-free computation of minmod of two numbers times 2s
static inline
float xmin2s(float s, float a, float b) {
    float sa = copysignf(s, a);
    float sb = copysignf(s, b);
    float abs_a = fabsf(a);
    float abs_b = fabsf(b);
    float min_abs = (abs_a < abs_b ? abs_a : abs_b);
    return (sa+sb) * min_abs;
}


// Limited combined slope estimate
static inline
float limdiff(float um, float u0, float up) {
    const float theta = 2.0;
    const float quarter = 0.25;
    float du1 = u0-um;   // Difference to left
    float du2 = up-u0;   // Difference to right
    float duc = up-um;   // Twice centered difference
    return xmin2s( quarter, xmin2s(theta, du1, du2), duc );
}


// Compute limited derivs
static inline
void limited_deriv1(float* restrict du,
                    const float* restrict u,
                    int ncell)
{
    for (int i = 0; i < ncell; ++i)
        du[i] = limdiff(u[i-1], u[i], u[i+1]);
}


// Compute limited derivs across stride
static inline
void limited_derivk(float* restrict du,
                    const float* restrict u,
                    int ncell, int stride)
{
    assert(stride > 0);
    for (int i = 0; i < ncell; ++i)
        du[i] = limdiff(u[i-stride], u[i], u[i+stride]);
}


/**
 * ### Advancing a time step
 *
 * Take one step of the numerical scheme.  This consists of two pieces:
 * a first-order corrector computed at a half time step, which is used
 * to obtain new $F$ and $G$ values; and a corrector step that computes
 * the solution at the full step.  For full details, we refer to the
 * [Jiang and Tadmor paper][jt].
 *
 * The `compute_step` function takes two arguments: the `io` flag
 * which is the time step modulo 2 (0 if even, 1 if odd); and the `dt`
 * flag, which actually determines the time step length.  We need
 * to know the even-vs-odd distinction because the Jiang-Tadmor
 * scheme alternates between a primary grid (on even steps) and a
 * staggered grid (on odd steps).  This means that the data at $(i,j)$
 * in an even step and the data at $(i,j)$ in an odd step represent
 * values at different locations in space, offset by half a space step
 * in each direction.  Every other step, we shift things back by one
 * mesh cell in each direction, essentially resetting to the primary
 * indexing scheme.
 *
 * We're slightly tricky in the corrector in that we write
 * $$
 *   v(i,j) = (s(i+1,j) + s(i,j)) - (d(i+1,j)-d(i,j))
 * $$
 * where $s(i,j)$ comprises the $u$ and $x$-derivative terms in the
 * update formula, and $d(i,j)$ the $y$-derivative terms.  This cuts
 * the arithmetic cost a little (not that it's that big to start).
 * It also makes it more obvious that we only need four rows worth
 * of scratch space.
 */


// Predictor half-step
static
void central2d_predict(float* restrict v,
                       float* restrict scratch,
                       const float* restrict u,
                       const float* restrict f,
                       const float* restrict g,
                       float dtcdx2, float dtcdy2,
                       int nx, int ny, int nfield)
{
    float* restrict fx = scratch;
    float* restrict gy = scratch+nx;
    for (int k = 0; k < nfield; ++k) {
        for (int iy = 1; iy < ny-1; ++iy) {
            int offset = (k*ny+iy)*nx+1;
            limited_deriv1(fx+1, f+offset, nx-2);
            limited_derivk(gy+1, g+offset, nx-2, nx);
            for (int ix = 1; ix < nx-1; ++ix) {
                int offset = (k*ny+iy)*nx+ix;
                v[offset] = u[offset] - dtcdx2 * fx[ix] - dtcdy2 * gy[ix];
            }
        }
    }
}


// Corrector
static
void central2d_correct_sd(float* restrict s,
                          float* restrict d,
                          const float* restrict ux,
                          const float* restrict uy,
                          const float* restrict u,
                          const float* restrict f,
                          const float* restrict g,
                          float dtcdx2, float dtcdy2,
                          int xlo, int xhi)
{
    for (int ix = xlo; ix < xhi; ++ix)
        s[ix] =
            0.2500f * (u [ix] + u [ix+1]) +
            0.0625f * (ux[ix] - ux[ix+1]) +
            dtcdx2  * (f [ix] - f [ix+1]);
    for (int ix = xlo; ix < xhi; ++ix)
        d[ix] =
            0.0625f * (uy[ix] + uy[ix+1]) +
            dtcdy2  * (g [ix] + g [ix+1]);
}


// Corrector
static
void central2d_correct(float* restrict v,
                       float* restrict scratch,
                       const float* restrict u,
                       const float* restrict f,
                       const float* restrict g,
                       float dtcdx2, float dtcdy2,
                       int xlo, int xhi, int ylo, int yhi,
                       int nx, int ny, int nfield)
{
    assert(0 <= xlo && xlo < xhi && xhi <= nx);
    assert(0 <= ylo && ylo < yhi && yhi <= ny);

    float* restrict ux = scratch;
    float* restrict uy = scratch +   nx;
    float* restrict s0 = scratch + 2*nx;
    float* restrict d0 = scratch + 3*nx;
    float* restrict s1 = scratch + 4*nx;
    float* restrict d1 = scratch + 5*nx;

    for (int k = 0; k < nfield; ++k) {

        float*       restrict vk = v + k*ny*nx;
        const float* restrict uk = u + k*ny*nx;
        const float* restrict fk = f + k*ny*nx;
        const float* restrict gk = g + k*ny*nx;

        limited_deriv1(ux+1, uk+ylo*nx+1, nx-2);
        limited_derivk(uy+1, uk+ylo*nx+1, nx-2, nx);
        central2d_correct_sd(s1, d1, ux, uy,
                             uk + ylo*nx, fk + ylo*nx, gk + ylo*nx,
                             dtcdx2, dtcdy2, xlo, xhi);

        for (int iy = ylo; iy < yhi; ++iy) {

            float* tmp;
            tmp = s0; s0 = s1; s1 = tmp;
            tmp = d0; d0 = d1; d1 = tmp;

            limited_deriv1(ux+1, uk+(iy+1)*nx+1, nx-2);
            limited_derivk(uy+1, uk+(iy+1)*nx+1, nx-2, nx);
            central2d_correct_sd(s1, d1, ux, uy,
                                 uk + (iy+1)*nx, fk + (iy+1)*nx, gk + (iy+1)*nx,
                                 dtcdx2, dtcdy2, xlo, xhi);

            for (int ix = xlo; ix < xhi; ++ix)
                vk[iy*nx+ix] = (s1[ix]+s0[ix])-(d1[ix]-d0[ix]);
        }
    }
}


static
void central2d_step(float* restrict u, float* restrict v,
                    float* restrict scratch,
                    float* restrict f,
                    float* restrict g,
                    int io, int nx, int ny, int ng,
                    int nfield, flux_t flux, speed_t speed,
                    float dt, float dx, float dy)
{
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;

    float dtcdx2 = 0.5 * dt / dx;
    float dtcdy2 = 0.5 * dt / dy;

    flux(f, g, u, nx_all * ny_all, nx_all * ny_all);

    central2d_predict(v, scratch, u, f, g, dtcdx2, dtcdy2,
                      nx_all, ny_all, nfield);

    // Flux values of f and g at half step
    for (int iy = 1; iy < ny_all-1; ++iy) {
        int jj = iy*nx_all+1;
        flux(f+jj, g+jj, v+jj, nx_all-2, nx_all * ny_all);
    }

    central2d_correct(v+io*(nx_all+1), scratch, u, f, g, dtcdx2, dtcdy2,
                      ng-io, nx+ng-io,
                      ng-io, ny+ng-io,
                      nx_all, ny_all, nfield);
}


/**
 * ### Advance a fixed time
 *
 * The `run` method advances from time 0 (initial conditions) to time
 * `tfinal`.  Note that `run` can be called repeatedly; for example,
 * we might want to advance for a period of time, write out a picture,
 * advance more, and write another picture.  In this sense, `tfinal`
 * should be interpreted as an offset from the time represented by
 * the simulator at the start of the call, rather than as an absolute time.
 *
 * We always take an even number of steps so that the solution
 * at the end lives on the main grid instead of the staggered grid.
 */

static
int central2d_xrun(central2d_t* restrict sim_, central2d_t* restrict sim, float tfinal)
{
   // Values from processor simulation domain
    float* restrict u = sim_ -> u;
    float* restrict v = sim_ -> v;
    float* restrict scratch = sim_ -> scratch;
    float* restrict f = sim_ -> f;
    float* restrict g = sim_ -> g;
    int nx = sim_ -> nx;
    int ny = sim_ -> ny;
    int ng = sim_ -> ng;
    int nfield = sim_ -> nfield;
    flux_t flux = sim_ -> flux;
    speed_t speed = sim_ -> speed;
    float dx = sim_ -> dx;
    float dy = sim_ -> dy;
    float cfl = sim_ -> cfl;
  
    int nstep = 0;
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
  
    bool done = false;
    float t = 0;
      printf("here1 \n");
      #pragma omp barrier
    while (!done) {
            printf("here2.5 \n");
          printf("here3 \n");
        float cxy[2] = {1.0e-15f, 1.0e-15f};
        speed(cxy, u, nx_all * ny_all, nx_all * ny_all);
        float dt = cfl / fmaxf(cxy[0]/dx, cxy[1]/dy);
                printf("here4 %d %d %d\n", nx, ny ,ng);
              printf("%f", u[5]);
          printf("blah");
        central2d_pushBC(u, sim->u,nx,ny,ng,sim->nx, sim->ny,sim->imin_, sim->jmin_); //push just BCs to global grid
        #pragma omp barrier //wait to make sure get most restrictive CFL
                      printf("here5 \n");
        central2d_BCset(u, sim->u, nx, ny, ng, sim->nx, sim->ny, sim_->imin_,sim_->jmin_,sim_->imax_,sim_->jmax_,
                        sim_->imin, sim->jmin, sim->imax, sim->jmax); //fill ghost cells from updated global grid
                            printf("here6 \n");
        if (t + 2*dt >= tfinal) {
            dt = (tfinal-t)/2;
            done = true;
        }
        central2d_step(u, v, scratch, f, g,
                       0, nx+4, ny+4, ng-2,
                       nfield, flux, speed,
                       dt, dx, dy);
        central2d_step(v, u, scratch, f, g,
                       1, nx, ny, ng,
                       nfield, flux, speed,
                       dt, dx, dy);
        t += 2*dt;
        nstep += 2;
    }
    return nstep;
}


int central2d_run(central2d_t* sim_, central2d_t* sim, float tfinal)
{
    return central2d_xrun(sim_,sim, tfinal);
}
