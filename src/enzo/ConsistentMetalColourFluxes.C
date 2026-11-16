/***********************************************************************
/
/  CONSISTENT METAL COLOUR FLUXES (species-tracked dust)
/
/  written by: Harrison Lo
/  date:       October 2026
/  modified1:  November 2026. The six components are the record and
/              Metal_Density is their sum; fraction reconstruction with
/              contact steepening replaces the face rescaling and the
/              outflow caps.
/
/  PURPOSE: With UseDustSpeciesTrack the gas-phase metals are carried as
/    six component colours (C, O, Mg, Si, Fe and the other metals) next to
/    Metal_Density. PPM reconstructs every colour on its own, so their
/    fluxes through a face do not add up, and the components drift away
/    from the total cell by cell. Here the six components are the record
/    and Metal_Density is kept equal to their sum. The x, y and z Euler
/    sweeps each do that in three steps:
/
/    1. ConsistentMetalDensity, before the reconstruction and again after
/       the Euler update: Metal_Density is set to the sum of the six on
/       the whole slice, ghost zones included. The reconstruction then
/       sees the total the components add up to, whatever interpolation,
/       flux correction or feedback did to the fields in between.
/
/    2. ConsistentMetalFaceStates, after inteuler: the six face states are
/       rebuilt from the mass fractions C/M. Each fraction gets a PPM
/       profile with intvar's slopes, monotonicity and flattening, and
/       contact steepening where the fraction jumps sharply but the gas
/       has no contact of its own (consistent multi-fluid advection,
/       Plewa & Mueller 1999, A&A 342, 179, section 2.3). The profiles are
/       traced like inteuler's colours, normalised to sum to one and
/       multiplied by the metal face state inteuler made, so the six face
/       states sum to it. Fractions keep a composition front two cells
/       wide where separately limited densities smear it over five.
/
/       Only the steepening of the Plewa & Mueller scheme is used. Its two
/       flattening steps (next to an extremum, and the coupled one that
/       makes the fractions sum to one before the trace) cost accuracy on
/       smooth profiles in 1D advection tests and gained nothing on fronts.
/
/    3. ConsistentMetalColourFluxes, after the Riemann solve: the metal
/       flux becomes the sum of the six component fluxes, so the Euler
/       update and the flux correction move exactly the metal that the
/       components move.
/
/    grid::SpeciesConsistent repeats step 1 on the whole grid before every
/    Grackle call. The counters below go to SpeciesConsistent.out with its
/    report.
/
/    Slice layout (all three sweeps): colour n occupies block n of
/    idim*jdim entries, element [j*idim + i], i along the sweep. Faces are
/    computed for i in [istart, iend+1] (0-based: GridStartIndex to
/    GridEndIndex+1); face i lies between cells i-1 and i.
/
************************************************************************/

#include <vector>
#include <cmath>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* Since the last report: 0 net and 1 absolute change of Metal_Density made
   by ConsistentMetalDensity on active cells (density x cell volume);
   2 fraction profiles built; 3 profiles steepened; 4 negative face
   fractions; 5 stencil cells without metal. */
static double CMFAccum[6] = {0, 0, 0, 0, 0, 0};

/* Colour numbers of Metal_Density and of the six components in the order
   C, O, Mg, Si, Fe, other. */

static int FindMetalColours(int NumberOfColours, int colnum[], int FieldType[],
                            int &nM, int six[6])
{
  const int types[6] = {MetalDensityCarbon, MetalDensityOxygen,
                        MetalDensityMagnesium, MetalDensitySilicon,
                        MetalDensityIron, MetalDensityOther};
  int n, q;
  nM = -1;
  for (q = 0; q < 6; q++)
    six[q] = -1;
  for (n = 0; n < NumberOfColours; n++) {
    int t = FieldType[colnum[n]];
    if (t == Metallicity)
      nM = n;
    for (q = 0; q < 6; q++)
      if (t == types[q])
        six[q] = n;
  }
  if (nM < 0)
    return FALSE;
  for (q = 0; q < 6; q++)
    if (six[q] < 0)
      return FALSE;
  return TRUE;
}

/* Step 1: Metal_Density = C + O + Mg + Si + Fe + other on the whole slice.
   A negative or non-finite component is an error. The change is counted
   on active cells only: i in [i1, i2], j in [j1, j2] of an active slice. */

int ConsistentMetalDensity(float *colslice, int NumberOfColours, int colnum[],
                           int FieldType[], int idim, int jdim,
                           int i1, int i2, int j1, int j2, int active,
                           double volume)
{
  if (!UseDustSpeciesTrack || NumberOfColours <= 0)
    return SUCCESS;
  int nM, six[6], q, i, j;
  if (!FindMetalColours(NumberOfColours, colnum, FieldType, nM, six))
    return FAIL;
  int size = idim*jdim;

  for (j = 0; j < jdim; j++)
    for (i = 0; i < idim; i++) {
      int m = j*idim + i;
      double sum = 0, old = colslice[nM*size + m];
      if (!std::isfinite(old))
        return FAIL;
      for (q = 0; q < 6; q++) {
        double v = colslice[six[q]*size + m];
        if (!std::isfinite(v) || v < 0)
          return FAIL;
        sum += v;
      }
      if (!std::isfinite(sum))
        return FAIL;
      colslice[nM*size + m] = (float) sum;
      if (active && i >= i1 && i <= i2 && j >= j1 && j <= j2) {
        double change = ((double) colslice[nM*size + m] - old) * volume;
        CMFAccum[0] += change;
        CMFAccum[1] += fabs(change);
      }
    }
  return SUCCESS;
}

static double Sign(double x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }
static double Clip(double x) { return max(0.0, min(1.0, x)); }

/* intvar's interpolation coefficients for cells a-2..b+2 on a possibly
   non-uniform dx: c[0], c[1] weight the centred slope and c[2]..c[5] build
   the edge value from the two neighbouring means and slopes. */

static void ProfileCoefficients(const float *dx, int a, int b,
                                std::vector<double> c[6])
{
  for (int i = a-2; i <= b+2; i++) {
    double t = dx[i]/(dx[i-1] + dx[i] + dx[i+1]);
    c[0][i] = t*(2*dx[i-1] + dx[i])/(dx[i+1] + dx[i]);
    c[1][i] = t*(2*dx[i+1] + dx[i])/(dx[i-1] + dx[i]);
  }
  for (int i = a-1; i <= b+2; i++) {
    double t = dx[i-2] + dx[i-1] + dx[i] + dx[i+1];
    double w = dx[i-1]/(dx[i-1] + dx[i]);
    double l = (dx[i-2] + dx[i-1])/(2*dx[i-1] + dx[i]);
    double r = (dx[i+1] + dx[i])/(2*dx[i] + dx[i-1]);
    w += 2*dx[i]*w/t*(l - r);
    c[2][i] = 1 - w;
    c[3][i] = w;
    c[4][i] = dx[i]/t*r;
    c[5][i] = -dx[i-1]/t*l;
  }
}

static void SecondDerivative(const double *q, const float *dx, int a, int b,
                             double *d2)
{
  for (int i = a-2; i <= b+2; i++)
    d2[i] = ((q[i+1] - q[i])/(dx[i+1] + dx[i]) -
             (q[i] - q[i-1])/(dx[i] + dx[i-1]))/(dx[i-1] + dx[i] + dx[i+1]);
}

/* PPM's contact-steepening strength in cell i (Colella & Woodward 1984,
   eq. 1.17, as in intvar): zero unless the second derivative changes sign
   across the cell and the jump exceeds 1% of the smaller neighbour. */

static double SteepeningStrength(const double *q, const float *dx,
                                 const double *d2, int i)
{
  double jump = q[i+1] - q[i-1];
  if (jump == 0 || d2[i+1]*d2[i-1] > 0 ||
      fabs(jump) <= .01*min(fabs(q[i+1]), fabs(q[i-1])))
    return 0;
  double dl = .5*(dx[i-1] + dx[i]), dr = .5*(dx[i] + dx[i+1]);
  double s1 = (d2[i-1] - d2[i+1])*(dl*dl*dl + dr*dr*dr)/((dl + dr)*jump);
  return Clip(20*(s1 - .05));
}

/* PPM profiles of the six fractions in one row. q holds the six rows of n
   fractions; the stencil needs q, rho and p in a-3..b+3 and the profiles
   (left and right edge values, and the curvature term) are returned for
   the donor cells a-1..b+1. */

static void FractionProfiles(const double *q, const float *rho, const float *p,
                             const float *dx, const float *flatten,
                             int n, int a, int b, double gamma, int iflatten,
                             const std::vector<double> c[6],
                             double *left, double *right, double *curve)
{
  std::vector<double> slope(n, 0), d2(n, 0), rg(n, 0), gd2(n, 0);
  std::vector<int> contact(n, 0);

  /* Cells where the gas itself has a contact discontinuity: a steepened
     density jump that is large next to the pressure jump. The gas
     reconstruction steepens those, so the fractions are left alone. */

  for (int i = a-3; i <= b+3; i++)
    rg[i] = rho[i];
  SecondDerivative(&rg[0], dx, a, b, &gd2[0]);
  for (int i = a-1; i <= b+1; i++) {
    double density_jump = fabs(rho[i+1] - rho[i-1])/min(rho[i+1], rho[i-1]);
    double pmin = min(p[i+1], p[i-1]);
    double pressure_jump = pmin > 0 ? fabs(p[i+1] - p[i-1])/pmin : HUGE_VAL;
    contact[i] = (SteepeningStrength(&rg[0], dx, &gd2[0], i) > 0 &&
                  gamma*.1*density_jump >= pressure_jump);
  }

  for (int s = 0; s < 6; s++) {
    const double *y = q + s*n;
    double *l = left + s*n, *r = right + s*n;

    /* Monotonised (van Leer) slopes and the edge values between cells. */

    for (int i = a-2; i <= b+2; i++) {
      double plus = y[i+1] - y[i], minus = y[i] - y[i-1];
      if (plus*minus > 0) {
        double cent = c[0][i]*plus + c[1][i]*minus;
        double vanl = 2*plus*minus/(minus + plus);
        slope[i] = Sign(cent)*min(min(fabs(cent), fabs(vanl)),
                                  min(2*fabs(minus), 2*fabs(plus)));
      } else
        slope[i] = 0;
    }
    for (int i = a-1; i <= b+2; i++) {
      l[i] = c[2][i]*y[i-1] + c[3][i]*y[i] + c[4][i]*slope[i-1] + c[5][i]*slope[i];
      r[i-1] = l[i];
    }

    /* Contact steepening: where the fraction jumps by more than three
       quarters of the change over the five-cell stencil, keeps its sense
       on both sides and the gas has no contact, move the edge values
       towards those extrapolated from the neighbours. */

    SecondDerivative(y, dx, a, b, &d2[0]);
    for (int i = a-1; i <= b+1; i++) {
      double den = y[i+2] - y[i-2];
      double eta = 0;
      if (den != 0 && (y[i+1] - y[i-1])/den > .75 &&
          (y[i+2] - y[i+1])*(y[i-1] - y[i-2]) > 0 && !contact[i])
        eta = SteepeningStrength(y, dx, &d2[0], i);
      if (eta > 0) {
        l[i] = (1 - eta)*l[i] + eta*(y[i-1] + .5*slope[i-1]);
        r[i] = (1 - eta)*r[i] + eta*(y[i+1] - .5*slope[i+1]);
        CMFAccum[3]++;
      }
    }

    /* intvar's order: monotonicity, shock flattening, then edge values
       bounded by the neighbouring means. */

    for (int i = a-1; i <= b+1; i++) {
      double t1 = (r[i] - y[i])*(y[i] - l[i]);
      double delta = r[i] - l[i], q6 = 6*(y[i] - .5*(r[i] + l[i]));
      if (t1 <= 0)
        l[i] = r[i] = y[i];
      if (delta*delta < delta*q6)
        l[i] = 3*y[i] - 2*r[i];
      if (delta*delta < -delta*q6)
        r[i] = 3*y[i] - 2*l[i];
      if (iflatten) {
        l[i] = y[i]*flatten[i] + l[i]*(1 - flatten[i]);
        r[i] = y[i]*flatten[i] + r[i]*(1 - flatten[i]);
      }
      l[i] = max(min(y[i], y[i-1]), min(max(y[i], y[i-1]), l[i]));
      r[i] = max(min(y[i], y[i+1]), min(max(y[i], y[i+1]), r[i]));
      CMFAccum[2]++;
    }
  }

  for (int i = a-1; i <= b+1; i++)
    for (int s = 0; s < 6; s++) {
      int k = s*n + i;
      curve[k] = 6*(q[k] - .5*(left[k] + right[k]));
    }
}

/* Step 2: face states of the six from their reconstructed fractions,
   normalised to the metal face state. js..je are the 0-based rows the
   sweep computes. Returns FAIL on invalid input (non-finite or negative
   values, non-positive gas density or pressure) and where the traced
   fractions do not sum to a positive number; a single negative traced
   fraction is passed on, and counted. */

int ConsistentMetalFaceStates(float *colls, float *colrs,
                              const float *colslice, const float *density,
                              const float *pressure, const float *velocity,
                              const float *dx, const float *flatten,
                              int NumberOfColours, int colnum[], int FieldType[],
                              int idim, int jdim, int istart, int iend,
                              int jstart, int jend, float dt, float gamma,
                              int pressure_free, int ordinary_flattening)
{
  if (!UseDustSpeciesTrack || NumberOfColours <= 0)
    return SUCCESS;
  int nM, six[6];
  if (!FindMetalColours(NumberOfColours, colnum, FieldType, nM, six))
    return FAIL;
  if (istart < 3 || iend+3 >= idim || iend < istart || jstart < 0 || jend >= jdim)
    return FAIL;
  for (int i = istart-3; i <= iend+3; i++)
    if (!(dx[i] > 0) || !std::isfinite(dx[i]))
      return FAIL;

  int size = idim*jdim;
  std::vector<double> q(6*idim, 0), l(6*idim, 0), r(6*idim, 0), q6(6*idim, 0);
  std::vector<double> coef[6];
  for (int k = 0; k < 6; k++)
    coef[k].resize(idim, 0);
  ProfileCoefficients(dx, istart, iend, coef);

  /* Sound speed inteuler uses for pressure-free flow. */
#ifdef CONFIG_PFLOAT_16
  const double pressure_free_sound_speed = 1.e-35;
#else
  const double pressure_free_sound_speed = 1.e-20;
#endif

  for (int j = jstart; j <= jend; j++) {
    int row = j*idim;

    /* Fractions C/M on the stencil. A cell without metal counts as all
       other metals, which keeps the fractions summing to one. */

    for (int i = istart-3; i <= iend+3; i++) {
      double M = colslice[nM*size + row + i], sum = 0;
      bool valid = std::isfinite(M) && M >= 0 && density[row+i] > 0 &&
        (pressure_free ? pressure[row+i] >= 0 : pressure[row+i] > 0) &&
        std::isfinite(density[row+i]) && std::isfinite(pressure[row+i]) &&
        std::isfinite(velocity[row+i]);
      for (int s = 0; s < 6; s++) {
        double C = colslice[six[s]*size + row + i];
        sum += C;
        valid = valid && std::isfinite(C) && C >= 0;
        q[s*idim + i] = M > 0 ? C/M : (s == 5 ? 1 : 0);
        valid = valid && std::isfinite(q[s*idim + i]);
      }
      if (M == 0) {
        CMFAccum[5]++;
        valid = valid && sum == 0;
      }
      if (!valid)
        return FAIL;
    }

    FractionProfiles(&q[0], density+row, pressure+row, dx, flatten+row,
                     idim, istart, iend, gamma, ordinary_flattening, coef,
                     &l[0], &r[0], &q6[0]);

    /* Trace each fraction to the face over the half step, as inteuler does
       for a colour: the left state of face i comes from cell i-1 and the
       right state from cell i, over the distance the fastest wave towards
       the face covers, or the flow itself when it moves towards the face. */

    for (int i = istart; i <= iend+1; i++)
      for (int side = 0; side < 2; side++) {
        int d = side == 0 ? i-1 : i;
        double u = velocity[row+d];
        double cs = pressure_free ? pressure_free_sound_speed :
          sqrt(gamma*pressure[row+d]/density[row+d]);
        double h = .5*dt/dx[d], c0 = h*u;
        double dist = side == 0 ? h*max(0.0, u + cs) : h*max(0.0, -u + cs);
        if ((side == 0 && u > 0) || (side == 1 && u < 0))
          dist = side == 0 ? c0 : -c0;
        double pred[6], sum = 0;
        for (int s = 0; s < 6; s++) {
          int k = s*idim + d;
          double delta = r[k] - l[k];
          pred[s] = side == 0 ? r[k] - dist*(delta - (1 - 4.0/3.0*dist)*q6[k]) :
                                l[k] + dist*(delta + (1 - 4.0/3.0*dist)*q6[k]);
          sum += pred[s];
          if (pred[s] < 0)
            CMFAccum[4]++;
        }
        float *out = side == 0 ? colls : colrs;
        double M = out[nM*size + row + i];
        if (!(sum > 0) || !std::isfinite(sum) || M < 0 || !std::isfinite(M))
          return FAIL;
        for (int s = 0; s < 6; s++)
          out[six[s]*size + row + i] = M*(pred[s]/sum);
      }
  }
  return SUCCESS;
}

/* Step 3: the metal flux is the sum of the six component fluxes. */

int ConsistentMetalColourFluxes(float *colf, int NumberOfColours, int colnum[],
                                int FieldType[], int idim, int jdim,
                                int istart, int iend)
{
  if (!UseDustSpeciesTrack || NumberOfColours <= 0)
    return SUCCESS;
  int nM, six[6], n, q, i, j;
  if (!FindMetalColours(NumberOfColours, colnum, FieldType, nM, six))
    return FAIL;
  int size = idim*jdim;

  /* Add the components in the order they are registered as colours. */
  int nsix = 0, order[6];
  for (n = 0; n < NumberOfColours; n++)
    for (q = 0; q < 6; q++)
      if (six[q] == n)
        order[nsix++] = n;

  for (j = 0; j < jdim; j++)
    for (i = istart; i <= iend+1; i++) {
      int m = j*idim + i;
      double sum = 0;
      for (q = 0; q < 6; q++)
        sum += colf[order[q]*size + m];
      colf[nM*size + m] = (float) sum;
    }
  return SUCCESS;
}

/* Hand the counters to SpeciesConsistentReport and reset them. */

void ConsistentMetalCounters(double values[6])
{
  for (int n = 0; n < 6; n++) {
    values[n] = CMFAccum[n];
    CMFAccum[n] = 0;
  }
}
