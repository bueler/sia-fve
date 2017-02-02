/* (C) 2015 Ed Bueler */

#include <stdarg.h>
#include <sys/time.h>
#include "io.h"

// set v = tol where v < tol
// (compare VecChop(), which sets v = 0 where abs(v) < tol)
PetscErrorCode VecTrueChop(Vec v, PetscReal tol) {
  PetscErrorCode  ierr;
  PetscReal       *a;
  PetscInt        n, i;

  ierr = VecGetLocalSize(v, &n); CHKERRQ(ierr);
  ierr = VecGetArray(v, &a); CHKERRQ(ierr);
  for (i = 0; i < n; ++i) {
      if (a[i] < tol)  a[i] = tol;
  }
  ierr = VecRestoreArray(v, &a); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// set v = 0 where w < tol
// v and w must have same layout
PetscErrorCode VecZeroWhereVecSmall(Vec v, Vec w, PetscReal tol) {
  PetscErrorCode  ierr;
  PetscReal       *av, *aw;
  PetscInt        n, i;

  ierr = VecGetLocalSize(v, &n); CHKERRQ(ierr);
  ierr = VecGetArray(v, &av); CHKERRQ(ierr);
  ierr = VecGetArray(w, &aw); CHKERRQ(ierr);
  for (i = 0; i < n; ++i) {
      if (aw[i] < tol)  av[i] = 0.0;
  }
  ierr = VecRestoreArray(v, &av); CHKERRQ(ierr);
  ierr = VecRestoreArray(w, &aw); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode AxisExtractSizes(Vec a, PetscInt *N, PetscReal *delta, PetscReal *L) {
    PetscErrorCode ierr;
    PetscReal *aa, fulllength;
    ierr = VecGetSize(a,N); CHKERRQ(ierr);
    if (*N < 4) {  // 4 is somewhat arbitrary
        SETERRQ(PETSC_COMM_WORLD,1,"read axis Vec has size too small\n");
    }
    ierr = VecGetArray(a, &aa);CHKERRQ(ierr);
    *delta = PetscAbsReal(aa[1] - aa[0]);
    fulllength = PetscAbsReal(aa[*N-1] - aa[0]) + *delta;
    *L = fulllength / 2.0;
    ierr = VecRestoreArray(a, &aa);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode ReadDimensions(AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer viewer;
    Vec x, y;

    PetscFunctionBeginUser;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,user->readname,FILE_MODE_READ,&viewer); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_SELF,&x); CHKERRQ(ierr);
    ierr = VecSetType(x, VECSEQ); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x,"x-axis from file"); CHKERRQ(ierr);
    ierr = VecLoad(x,viewer); CHKERRQ(ierr);
    ierr = AxisExtractSizes(x,&user->Nx,&user->dx,&user->Lx); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_SELF,&y); CHKERRQ(ierr);
    ierr = VecSetType(y, VECSEQ); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)y,"y-axis from file"); CHKERRQ(ierr);
    ierr = VecLoad(y,viewer); CHKERRQ(ierr);
    ierr = AxisExtractSizes(y,&user->Ny,&user->dy,&user->Ly); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode ReadAndReshape2DVec(Vec v, PetscViewer viewer, AppCtx *user) {
    PetscErrorCode ierr;
    Vec            ser;
    PetscInt       j, k, s,
                   NN = user->Nx * user->Ny;
    PetscReal      *aser, **av;
    DMDALocalInfo  info;

    PetscFunctionBeginUser;
    ierr = VecCreateSeq(PETSC_COMM_SELF,NN,&ser); CHKERRQ(ierr);
    ierr = VecLoad(ser,viewer); CHKERRQ(ierr);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = VecGetArray(ser, &aser);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, v, &av);CHKERRQ(ierr);
    for (k=info.ys; k<info.ys+info.ym; k++) {
        for (j=info.xs; j<info.xs+info.xm; j++) {
            s = k * info.mx + j;  //  0 <= s <= info.my*info.mx - 1
            av[k][j] = aser[s];
        }
    }
    ierr = VecRestoreArray(ser, &aser);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, v, &av);CHKERRQ(ierr);

    ierr = VecDestroy(&ser); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode DiscardDimensions(PetscViewer viewer) {
    PetscErrorCode ierr;
    Vec            x, y;
    ierr = VecCreate(PETSC_COMM_SELF,&x); CHKERRQ(ierr);
    ierr = VecLoad(x,viewer); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_SELF,&y); CHKERRQ(ierr);
    ierr = VecLoad(y,viewer); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode DiscardVecs(PetscInt n, PetscViewer viewer, AppCtx *user) {
    PetscErrorCode ierr;
    Vec            tmp;
    PetscInt       j;
    for (j = 0; j < n; j++) {
        ierr = VecDuplicate(user->b,&tmp); CHKERRQ(ierr);
        ierr = ReadAndReshape2DVec(tmp, viewer, user); CHKERRQ(ierr);
        ierr = VecDestroy(&tmp); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode ReadDataVecs(AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer    viewer;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,user->readname,FILE_MODE_READ,&viewer); CHKERRQ(ierr);
    ierr = DiscardDimensions(viewer); CHKERRQ(ierr);

    ierr = ReadAndReshape2DVec(user->b, viewer, user); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(user->m, viewer, user); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(user->Hexact, viewer, user); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(user->Hinitial, viewer, user); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode ReadInitialH(AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer    viewer;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,user->readinitialname,FILE_MODE_READ,&viewer); CHKERRQ(ierr);
    ierr = DiscardDimensions(viewer); CHKERRQ(ierr);
    // read and discard b,m,Hexact
    ierr = DiscardVecs(3,viewer,user); CHKERRQ(ierr);
    // actually read Hinitial
    ierr = ReadAndReshape2DVec(user->Hinitial, viewer, user); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode GenerateInitialHFromReadSurface(AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer    viewer;
    Vec            tmpb, tmpH;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,user->readinitialname,FILE_MODE_READ,&viewer); CHKERRQ(ierr);
    ierr = DiscardDimensions(viewer); CHKERRQ(ierr);
    // read b into tmpb
    ierr = VecDuplicate(user->b,&tmpb); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(tmpb, viewer, user); CHKERRQ(ierr);
    // read and discard m and Hexact
    ierr = DiscardVecs(2,viewer,user); CHKERRQ(ierr);
    // read H into tmpH
    ierr = VecDuplicate(user->Hinitial,&tmpH); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(tmpH, viewer, user); CHKERRQ(ierr);
    // add tmpb so Hinitial is really s
    ierr = VecWAXPY(user->Hinitial,1.0,tmpb,tmpH); CHKERRQ(ierr);  // Hinitial <- 1.0*tmpb + tmpH
    // subtract previously-read b to create new thickness, but note this is
    //   nonsense in ice-free areas; it is just the difference of beds;
    //   and it can come out negative even where tmpH > 0
    ierr = VecAXPY(user->Hinitial,-1.0,user->b); CHKERRQ(ierr);  // Hinitial <- -1.0*b + Hinitial
    // now go though Hinitial and zero it out everywhere read tmpH is nearly zero (tmpH<1.0)
    ierr = VecZeroWhereVecSmall(user->Hinitial,tmpH,1.0); CHKERRQ(ierr);
    // zero out negative values
    ierr = VecTrueChop(user->Hinitial,0.0); CHKERRQ(ierr);
    ierr = VecDestroy(&tmpb); CHKERRQ(ierr);
    ierr = VecDestroy(&tmpH); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode DumpToFile(Vec H, Vec r, const char name[], AppCtx *user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            x, y;
    PetscInt       j, k;
    PetscReal      *ax, *ay;
    char           filename[1024];
    int            strerr;
    PetscViewer    viewer;

    if (!user->silent) PetscPrintf(PETSC_COMM_WORLD,"writing x,y,b,m,Hexact,H,residual%s into file %s in %s ...\n",
             (user->nodiag) ? "" : ",D,Wmag",name,user->figsprefix);

    strerr = sprintf(filename,"%s%s",user->figsprefix,name);
    if (strerr < 0) {
        SETERRQ1(PETSC_COMM_WORLD,6,"sprintf() returned %d < 0 ... stopping\n",strerr);
    }

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,info.mx,&x);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,info.my,&y);CHKERRQ(ierr);
    ierr = VecGetArray(x, &ax);CHKERRQ(ierr);
    for (j=0; j<info.mx; j++) {
        ax[j] = -user->Lx + user->dx/2.0 + j * user->dx;
    }
    ierr = VecRestoreArray(x, &ax);CHKERRQ(ierr);
    ierr = VecGetArray(y, &ay);CHKERRQ(ierr);
    for (k=0; k<info.my; k++) {
        ay[k] = -user->Ly + user->dy/2.0 + k * user->dy;
    }
    ierr = VecRestoreArray(y, &ay);CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
    ierr = VecView(x,viewer); CHKERRQ(ierr);
    ierr = VecView(y,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_APPEND,&viewer); CHKERRQ(ierr);
    ierr = VecView(user->b,viewer); CHKERRQ(ierr);
    ierr = VecView(user->m,viewer); CHKERRQ(ierr);
    ierr = VecView(user->Hexact,viewer); CHKERRQ(ierr);
    ierr = VecView(H,viewer); CHKERRQ(ierr);

    // keep writing more than
    ierr = VecView(r,viewer); CHKERRQ(ierr);
    if (!user->nodiag) {
        DiagnosticScheme *ds = &(user->ds);
        ierr = VecView(ds->Dnodemax,viewer); CHKERRQ(ierr);
        ierr = VecView(ds->Wmagnodemax,viewer); CHKERRQ(ierr);
    }

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode GetVolumeArea(Vec H, AppCtx *user, PetscReal *volH, PetscReal *volHexact, PetscReal *areaH) {
  PetscErrorCode  ierr;
  const PetscReal darea = user->dx * user->dy;
  PetscReal       *aH;
  PetscInt        n, i, loccount=0, count;

  ierr = VecSum(H,volH); CHKERRQ(ierr);
  *volH *= darea;

  ierr = VecSum(user->Hexact,volHexact); CHKERRQ(ierr);
  *volHexact *= darea;

  VecGetLocalSize(H, &n);
  VecGetArray(H, &aH);
  for (i = 0; i < n; ++i) {
      if (aH[i] > 1.0)   // 1 m tolerance; fails if "> 0.0"
          loccount += 1;
  }
  VecRestoreArray(H, &aH);
  ierr = MPI_Allreduce(&loccount,&count,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
  *areaH = (PetscReal)count * darea;

  PetscFunctionReturn(0);
}

PetscErrorCode GetErrors(Vec H, AppCtx *user, PetscReal *enorminf, PetscReal *enorm1) {
  PetscErrorCode  ierr;
  Vec             dH;
  ierr = VecDuplicate(H,&dH); CHKERRQ(ierr);
  ierr = VecWAXPY(dH,-1.0,user->Hexact,H); CHKERRQ(ierr);  // dH := (-1.0) Hexact + H = H - Hexact
  ierr = VecNorm(dH,NORM_INFINITY,enorminf); CHKERRQ(ierr);
  ierr = VecNorm(dH,NORM_1,enorm1); CHKERRQ(ierr);
  ierr = VecDestroy(&dH); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DiffusivityReduce(AppCtx *user, PetscReal *avD, PetscReal *maxD) {
  PetscErrorCode   ierr;
  DiagnosticScheme *ds = &(user->ds);
  PetscInt         avDcount;
  ierr = MPI_Allreduce(&ds->avD,avD,1,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = MPI_Allreduce(&ds->avDcount,&avDcount,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
  *avD /= avDcount;
  ierr = MPI_Allreduce(&ds->maxD,maxD,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StdoutReport(Vec H, AppCtx *user) {
  PetscErrorCode  ierr;
  DMDALocalInfo   info;
  PetscInt        NN;
  PetscReal       volH, volHexact, areaH, enorminf, enorm1, voldiffrel;

  ierr = GetVolumeArea(H, user, &volH, &volHexact, &areaH); CHKERRQ(ierr);
  if (!user->silent) PetscPrintf(PETSC_COMM_WORLD,"       vol = %8.2e km3,  area = %8.2e km2",
                (double)volH / 1.0e9, (double)areaH / 1.0e6);

  if (!user->nodiag) {
      PetscReal avD, maxD;
      ierr = DiffusivityReduce(user,&avD,&maxD); CHKERRQ(ierr);
      if (!user->silent) PetscPrintf(PETSC_COMM_WORLD,";  diagnostics:  max D = %6.4f,  av D = %6.4f m^2 s-1\n",
                    (double)maxD, (double)avD);
  } else
      if (!user->silent) PetscPrintf(PETSC_COMM_WORLD,"\n");

  ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
  NN = info.mx * info.my;
  ierr = GetErrors(H, user, &enorminf, &enorm1); CHKERRQ(ierr);
  voldiffrel = PetscAbsReal(volH - volHexact) / volHexact;
  if (!user->silent) PetscPrintf(PETSC_COMM_WORLD,"       errors:  max = %7.2f m,  av = %7.2f m,  voldiff%% = %5.2f\n",
                (double)enorminf, (double)enorm1 / NN, 100.0 * (double)voldiffrel);

  PetscFunctionReturn(0);
}

PetscErrorCode WriteHistoryFile(Vec H, const char name[], int argc, char **argv, AppCtx *user) {
    PetscErrorCode  ierr;
    PetscViewer     viewer;
    DMDALocalInfo   info;
    PetscInt        NN;
    char            cmdline[1024] = "", filename[1024] = "";
    int             j, strerr, size;
    double          computationtime;
    PetscReal       volH, volHexact, areaH, enorminf, enorm1, avD, maxD;

    if (!user->silent) PetscPrintf(PETSC_COMM_WORLD,"writing %s to %s ...\n",name,user->figsprefix);

    strerr = sprintf(filename,"%s%s",user->figsprefix,name);
    if (strerr < 0) { SETERRQ1(PETSC_COMM_WORLD,6,"sprintf() returned %d < 0 ... stopping\n",strerr); }
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer); CHKERRQ(ierr);
    for (j=0; j<argc; j++) {
        strcat(cmdline,argv[j]);
        strcat(cmdline," ");
    }
    ierr = PetscViewerASCIIPrintf(viewer,"%s\n",cmdline); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
    NN = info.mx * info.my;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"number of processors  %d\n",size); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"grid points in x-direction  %d\n",info.mx); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"grid points in y-direction  %d\n",info.my); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"domain half-width in x-direction (m)  %.6f\n",(double)user->Lx); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"domain half-width in y-direction (m)  %.6f\n",(double)user->Ly); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"spacing in x-direction (m)  %.6f\n",(double)user->dx); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"spacing in y-direction (m)  %.6f\n",(double)user->dy); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"last successful value of eps  %.6e\n",(double)user->eps); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"last c.s. level where convergence happened  %d\n",user->goodm); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"count of recovery steps (zero if no recovery)  %d\n",user->recoverycount); CHKERRQ(ierr);
    ierr = GetVolumeArea(H, user, &volH, &volHexact, &areaH); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution ice volume (m^3)  %.6e\n",(double)volH); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"exact ice volume (m^3)  %.6e\n",(double)volHexact); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution ice area (m^2)  %.6e\n",(double)areaH); CHKERRQ(ierr);
    ierr = DiffusivityReduce(user,&avD,&maxD); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"average solution diffusivity (m^2 s-1)  %.6e\n",(double)avD); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"maximum solution diffusivity (m^2 s-1)  %.6e\n",(double)maxD); CHKERRQ(ierr);
    ierr = GetErrors(H, user, &enorminf, &enorm1); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"max thickness error (m)  %.6e\n",(double)enorminf); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"av thickness error (m)  %.6e\n",(double)enorm1 / NN); CHKERRQ(ierr);
    computationtime =   (double)(user->endtime.tv_usec - user->starttime.tv_usec) / 1.0e6
                      + (double)(user->endtime.tv_sec -  user->starttime.tv_sec);
    ierr = PetscViewerASCIIPrintf(viewer,"total time (seconds)  %.6f\n",computationtime); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode ShowOne(Vec v, PetscInt xdim, PetscInt ydim, const char *title) {
  PetscErrorCode ierr;
  PetscViewer    graphical;
  PetscFunctionBeginUser;
  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,title,
                             PETSC_DECIDE,PETSC_DECIDE,xdim,ydim,&graphical); CHKERRQ(ierr);
  ierr = VecView(v,graphical); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&graphical); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ShowFields(AppCtx *user) {
  PetscErrorCode ierr;
  PetscInt       xdim, ydim;
  DMDALocalInfo  info;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  xdim = PetscMax(200,PetscMin(300,info.mx));
  ydim = PetscMax(200,PetscMin(561,info.my));
  ierr = ShowOne(user->b,xdim,ydim,"bed topography (m)"); CHKERRQ(ierr);
  ierr = ShowOne(user->m,xdim,ydim,"climatic mass balance (m s-1)"); CHKERRQ(ierr);
  ierr = ShowOne(user->Hexact,xdim,ydim,"exact or observed thickness (m)"); CHKERRQ(ierr);
  ierr = ShowOne(user->Hinitial,xdim,ydim,"initial thickness (m)"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

