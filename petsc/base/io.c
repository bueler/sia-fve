/* (C) 2015 Ed Bueler */

#include <stdarg.h>
#include <sys/time.h>
#include "io.h"

void myPrintf(const AppCtx *user, const char format[], ...) {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if ((!user->silent) && (rank==0)) {
        va_list Argp;
        va_start(Argp, format);
        PetscVFPrintf(PETSC_STDOUT, format, Argp);
        va_end(Argp);
    }
}

PetscErrorCode ReadDimensions(AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer viewer;
    Vec x, y;
    PetscReal *ax, *ay, fulllengthx, fulllengthy;

    PetscFunctionBeginUser;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,user->readname,FILE_MODE_READ,&viewer); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_SELF,&x); CHKERRQ(ierr);
    ierr = VecSetType(x, VECSEQ); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x,"x-axis from file"); CHKERRQ(ierr);
    ierr = VecLoad(x,viewer); CHKERRQ(ierr);
    ierr = VecGetSize(x,&user->Nx); CHKERRQ(ierr);
    if (user->Nx < 4) {  // 4 is somewhat arbitrary
        SETERRQ(PETSC_COMM_WORLD,1,"read Vec x has size too small\n");
    }
    ierr = VecGetArray(x, &ax);CHKERRQ(ierr);
    user->dx = PetscAbsReal(ax[1] - ax[0]);
    fulllengthx = PetscAbsReal(ax[user->Nx-1] - ax[0]) + user->dx;
    user->Lx = fulllengthx / 2.0;
    ierr = VecRestoreArray(x, &ax);CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_SELF,&y); CHKERRQ(ierr);
    ierr = VecSetType(y, VECSEQ); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)y,"y-axis from file"); CHKERRQ(ierr);
    ierr = VecLoad(y,viewer); CHKERRQ(ierr);
    ierr = VecGetSize(y,&user->Ny); CHKERRQ(ierr);
    if (user->Ny < 4) {  // 4 is somewhat arbitrary
        SETERRQ(PETSC_COMM_WORLD,2,"read Vec y has size too small\n");
    }
    ierr = VecGetArray(y, &ay);CHKERRQ(ierr);
    fulllengthy = PetscAbsReal(ay[user->Ny-1] - ay[0]) + user->dx;
    user->Ly = fulllengthy / 2.0;
    ierr = VecRestoreArray(y, &ay);CHKERRQ(ierr);
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


// NOTE: order in which data is read must match order in which nc2petsc.py writes
PetscErrorCode ReadDataVecs(AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer    viewer;
    Vec            x, y;  // will be read and discarded

    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,user->readname,FILE_MODE_READ,&viewer); CHKERRQ(ierr);

    // read and discard x and y
    ierr = VecCreate(PETSC_COMM_SELF,&x); CHKERRQ(ierr);
    ierr = VecLoad(x,viewer); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_SELF,&y); CHKERRQ(ierr);
    ierr = VecLoad(y,viewer); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);

    ierr = ReadAndReshape2DVec(user->b, viewer, user); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(user->m, viewer, user); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(user->Hexact, viewer, user); CHKERRQ(ierr);
    ierr = ReadAndReshape2DVec(user->Hinitial, viewer, user); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


//  write a vector into a petsc binary file with name built from prefix and name
PetscErrorCode ViewToBinary(PetscBool fromrankzero, Vec v, const char prefix[], const char name[]) {
    PetscErrorCode ierr;
    PetscViewer    viewer;
    char           filename[1024];
    int            strerr;
    strerr = sprintf(filename,"%s%s",prefix,name);
    if (strerr < 0) {
        SETERRQ1(PETSC_COMM_WORLD,6,"sprintf() returned %d < 0 ... stopping\n",strerr);
    }
    if (fromrankzero == PETSC_TRUE) {
        int  rank;
        ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
        if (rank == 0) {
            ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
            ierr = VecView(v,viewer); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        }
    } else {
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
        ierr = VecView(v,viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}


//  write various vectors to petsc binary files, one file per vec
//  FIXME: this could be re-designed to write a combined file like the one
//         we read; this would require modifications to figsmahaffy.py
PetscErrorCode DumpToFiles(Vec H, Vec r, AppCtx *user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            x, y;
    PetscInt       j, k;
    PetscReal      *ax, *ay;

    myPrintf(user,"writing {x,y,H,b,m,Hexact,residual%s}.dat to %s ...\n",
             (user->nodiag) ? "" : ",D,Wmag",user->figsprefix);

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
        ay[k] = -user->Ly + user->dx/2.0 + k * user->dx;
    }
    ierr = VecRestoreArray(y, &ay);CHKERRQ(ierr);
    ierr = ViewToBinary(PETSC_TRUE,x,user->figsprefix,"x.dat"); CHKERRQ(ierr);
    ierr = ViewToBinary(PETSC_TRUE,y,user->figsprefix,"y.dat"); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);

    ierr = ViewToBinary(PETSC_FALSE,user->b,user->figsprefix,"b.dat"); CHKERRQ(ierr);
    ierr = ViewToBinary(PETSC_FALSE,user->m,user->figsprefix,"m.dat"); CHKERRQ(ierr);
    ierr = ViewToBinary(PETSC_FALSE,user->Hexact,user->figsprefix,"Hexact.dat"); CHKERRQ(ierr);
    ierr = ViewToBinary(PETSC_FALSE,H,user->figsprefix,"H.dat"); CHKERRQ(ierr);
    ierr = ViewToBinary(PETSC_FALSE,r,user->figsprefix,"residual.dat"); CHKERRQ(ierr);

    if (!user->nodiag) {
        DiagnosticScheme *ds = &(user->ds);
        ierr = ViewToBinary(PETSC_FALSE,ds->Dnodemax,user->figsprefix,"D.dat"); CHKERRQ(ierr);
        ierr = ViewToBinary(PETSC_FALSE,ds->Wmagnodemax,user->figsprefix,"Wmag.dat"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


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


// volH and volHexact have units m^3; areaH has units m^2
PetscErrorCode GetVolumeArea(Vec H, AppCtx *user, PetscReal *volH, PetscReal *volHexact, PetscReal *areaH) {
  PetscErrorCode  ierr;
  const PetscReal darea = user->dx * user->dx;
  PetscReal      *aH;
  PetscInt       n, i, loccount=0, count;

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


// enorminf = ||H - Hexact||_infty,  enorm1 = ||H-Hexact||_1;  both have units of m
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


PetscErrorCode StdoutReport(Vec H, AppCtx *user) {
  PetscErrorCode  ierr;
  DMDALocalInfo   info;
  PetscInt        NN;
  PetscReal       volH, volHexact, areaH, enorminf, enorm1, voldiffrel;

  ierr = GetVolumeArea(H, user, &volH, &volHexact, &areaH); CHKERRQ(ierr);
  myPrintf(user,"       vol = %8.2e km3,  area = %8.2e km2",
                (double)volH / 1.0e9, (double)areaH / 1.0e6);

  if (!user->nodiag) {
      DiagnosticScheme *ds = &(user->ds);
      PetscInt  avDcount;
      PetscReal avD, maxD;
      ierr = MPI_Allreduce(&ds->avD,&avD,1,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
      ierr = MPI_Allreduce(&ds->avDcount,&avDcount,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
      avD /= avDcount;
      ierr = MPI_Allreduce(&ds->maxD,&maxD,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);
      myPrintf(user,";  diagnostics:  max D = %6.4f,  av D = %6.4f m^2 s-1\n",
                    (double)maxD, (double)avD);
  } else
      myPrintf(user,"\n");

  ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
  NN = info.mx * info.my;
  ierr = GetErrors(H, user, &enorminf, &enorm1); CHKERRQ(ierr);
  voldiffrel = PetscAbsReal(volH - volHexact) / volHexact;
  myPrintf(user,"       errors:  max = %7.2f m,  av = %7.2f m,  voldiff%% = %5.2f\n",
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
    PetscReal       volH, volHexact, areaH, enorminf, enorm1;

    myPrintf(user,"writing %s to %s ...\n",name,user->figsprefix);

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
    ierr = PetscViewerASCIIPrintf(viewer,"spacing in y-direction (m)  %.6f\n",(double)user->dx); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"last successful value of eps  %.6e\n",(double)user->eps); CHKERRQ(ierr);
    ierr = GetVolumeArea(H, user, &volH, &volHexact, &areaH); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution ice volume (m^3)  %.6e\n",(double)volH); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"exact ice volume (m^3)  %.6e\n",(double)volHexact); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution ice area (m^2)  %.6e\n",(double)areaH); CHKERRQ(ierr);
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

