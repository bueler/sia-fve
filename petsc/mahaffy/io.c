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


PetscErrorCode WriteHistoryFile(const char name[],int argc,char **argv, AppCtx *user) {
    SETERRQ(PETSC_COMM_WORLD,1,"NOT IMPLEMENTED ... stopping\n");
/*
struct timeval  tv1, tv2;
gettimeofday(&tv1, NULL);
// stuff
gettimeofday(&tv2, NULL);
printf ("Total time = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
*/
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

