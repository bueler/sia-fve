/* (C) 2015 Ed Bueler */

#include "binaryio.h"

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
PetscErrorCode DumpToFiles(Vec H, AppCtx *user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            x, y, dH;
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
    ierr = ViewToBinary(PETSC_FALSE,H,user->figsprefix,"H.dat"); CHKERRQ(ierr);

    ierr = VecDuplicate(H,&dH); CHKERRQ(ierr);
    ierr = VecWAXPY(dH,-1.0,user->Hexact,H); CHKERRQ(ierr);    // dH := (-1.0) Hexact + H = H - Hexact
    ierr = ViewToBinary(PETSC_FALSE,dH,user->figsprefix,"Herror.dat"); CHKERRQ(ierr);
    ierr = VecDestroy(&dH); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


PetscErrorCode CreateReadVecs(SerialReadVecs *readvecs,AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer viewer;
    Vec x, y;
    PetscReal *ax, *ay, fulllengthx, fulllengthy;
    PetscInt  NN;

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

    NN = user->Nx * user->Ny;

    ierr = VecCreateSeq(PETSC_COMM_SELF,NN,&readvecs->topgread); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)readvecs->topgread,"bed topography from file (serial)"); CHKERRQ(ierr);
    ierr = VecLoad(readvecs->topgread,viewer); CHKERRQ(ierr);

    ierr = VecCreateSeq(PETSC_COMM_SELF,NN,&readvecs->cmbread); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)readvecs->cmbread,"climatic mass balance from file (serial)"); CHKERRQ(ierr);
    ierr = VecLoad(readvecs->cmbread,viewer); CHKERRQ(ierr);

    ierr = VecCreateSeq(PETSC_COMM_SELF,NN,&readvecs->thkobsread); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)readvecs->thkobsread,"observed thickness from file (serial)"); CHKERRQ(ierr);
    ierr = VecLoad(readvecs->thkobsread,viewer); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


PetscErrorCode ReshapeAndDestroyReadVecs(SerialReadVecs *readvecs, AppCtx *user) {
    PetscErrorCode ierr;
    PetscInt       j, k, s;
    PetscReal      *atopg, *acmb, *athkobs, **ab, **am, **aHexact;
    DMDALocalInfo  info;

    PetscFunctionBeginUser;
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = VecGetArray(readvecs->topgread, &atopg);CHKERRQ(ierr);
    ierr = VecGetArray(readvecs->cmbread, &acmb);CHKERRQ(ierr);
    ierr = VecGetArray(readvecs->thkobsread, &athkobs);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->b, &ab);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->m, &am);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->Hexact, &aHexact);CHKERRQ(ierr);
    for (k=info.ys; k<info.ys+info.ym; k++) {
        for (j=info.xs; j<info.xs+info.xm; j++) {
            s = k * info.mx + j;  //  0 <= s <= info.my*info.mx - 1
            ab[k][j] = atopg[s];
            am[k][j] = acmb[s];
            aHexact[k][j] = athkobs[s];
        }
    }
    ierr = VecRestoreArray(readvecs->topgread, &atopg);CHKERRQ(ierr);
    ierr = VecRestoreArray(readvecs->cmbread, &acmb);CHKERRQ(ierr);
    ierr = VecRestoreArray(readvecs->thkobsread, &athkobs);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->b, &ab);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->m, &am);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->Hexact, &aHexact);CHKERRQ(ierr);

    ierr = VecDestroy(&readvecs->topgread); CHKERRQ(ierr);
    ierr = VecDestroy(&readvecs->cmbread); CHKERRQ(ierr);
    ierr = VecDestroy(&readvecs->thkobsread); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

