/* (C) 2015 Ed Bueler */

#include <sys/time.h>
#include "io.h"

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
    user->dx = ax[1] - ax[0];
    fulllengthx = ax[user->Nx-1] - ax[0] + user->dx;
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
    fulllengthy = ay[user->Ny-1] - ay[0] + user->dx;
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
//         we read, but with H and Herror appended, but this would require
//         mods to figsmahaffy.py
PetscErrorCode DumpToFiles(Vec H, AppCtx *user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            x, y;
    PetscInt       j, k;
    PetscReal      *ax, *ay;

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

    PetscFunctionReturn(0);
}


// volH and volHexact have units m^3
PetscErrorCode GetVolumes(Vec H, AppCtx *user, PetscReal *volH, PetscReal *volHexact) {
  PetscErrorCode  ierr;
  const PetscReal darea = user->dx * user->dx;
  ierr = VecSum(H,volH); CHKERRQ(ierr);
  ierr = VecSum(user->Hexact,volHexact); CHKERRQ(ierr);
  *volH      *= darea;
  *volHexact *= darea;
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
  PetscReal       maxD, volH, volHexact, enorminf, enorm1, voldiffrel;

  ierr = MPI_Allreduce(&user->maxD,&maxD,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = GetVolumes(H, user, &volH, &volHexact); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
             "        state:  vol = %8.4e km^3,  max D = %8.4f m^2 s-1\n",
             (double)volH / 1.0e9, (double)maxD); CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
  NN = info.mx * info.my;
  ierr = GetErrors(H, user, &enorminf, &enorm1); CHKERRQ(ierr);
  voldiffrel = PetscAbsReal(volH - volHexact) / volHexact;
  ierr = PetscPrintf(PETSC_COMM_WORLD,
             "       errors:  max = %7.2f m,       av = %7.2f m,        voldiff%% = %5.2f\n",
             (double)enorminf, (double)enorm1 / NN, 100.0 * (double)voldiffrel); CHKERRQ(ierr);
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
    PetscReal       volH, volHexact, enorminf, enorm1;

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
    ierr = PetscViewerASCIIPrintf(viewer,"last successful value of slopeeps  %.6e\n",(double)user->slopeeps); CHKERRQ(ierr);
    ierr = GetVolumes(H, user, &volH, &volHexact); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution ice volume (m^3)  %.6e\n",(double)volH); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"exact ice volume (m^3)  %.6e\n",(double)volHexact); CHKERRQ(ierr);
    ierr = GetErrors(H, user, &enorminf, &enorm1); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"max thickness error (m)  %.6e\n",(double)enorminf); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"av thickness error (m)  %.6e\n",(double)enorm1 / NN); CHKERRQ(ierr);
    computationtime =   (double)(user->endtime.tv_usec - user->starttime.tv_usec) / 1.0e6
                      + (double)(user->endtime.tv_sec -  user->starttime.tv_sec);
    ierr = PetscViewerASCIIPrintf(viewer,"total time (seconds)  %.6f\n",computationtime); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


PetscErrorCode ShowFields(AppCtx *user) {
  PetscErrorCode ierr;
  PetscViewer    graphical;
  PetscInt       xdim, ydim;
  DMDALocalInfo  info;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  xdim = PetscMax(200,PetscMin(300,info.mx));
  ydim = PetscMax(200,PetscMin(561,info.my));

  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"bed topography (m)",
                             PETSC_DECIDE,PETSC_DECIDE,xdim,ydim,&graphical); CHKERRQ(ierr);
  ierr = VecView(user->b,graphical); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&graphical); CHKERRQ(ierr);
  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"climatic mass balance (m s-1)",
                             PETSC_DECIDE,PETSC_DECIDE,xdim,ydim,&graphical); CHKERRQ(ierr);
  ierr = VecView(user->m,graphical); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&graphical); CHKERRQ(ierr);
  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"exact or observed thickness (m)",
                             PETSC_DECIDE,PETSC_DECIDE,xdim,ydim,&graphical); CHKERRQ(ierr);
  ierr = VecView(user->Hexact,graphical); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&graphical); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

