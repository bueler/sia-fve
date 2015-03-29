/* I don't want to lose this idea, but currently it has no use. */

PetscErrorCode ExplicitStepSmoother(Vec H, AppCtx *user) {
    PetscErrorCode  ierr;
    DMDALocalInfo   info;
    PetscReal       **aHloc, **aR, **aH, maxD, tmpeps, deltat, mu;
    PetscInt        j, k;
    Vec             Hloc, R;

    PetscFunctionBeginUser;
    // generate ghosted version of thickness H
    ierr = DMCreateLocalVector(user->da, &Hloc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da,H,INSERT_VALUES,Hloc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da,H,INSERT_VALUES,Hloc); CHKERRQ(ierr);

    // prepare to get residual from FormFunctionLocal()
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->da, &R); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, Hloc, &aHloc);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, H, &aH);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, R, &aR);CHKERRQ(ierr);

    // get eps=0 residual corresponding to H
    //user->eps = 0.0;
    tmpeps = user->eps; // set aside current value of eps
    ierr = FormFunctionLocal(&info,aHloc,aR,user);CHKERRQ(ierr); // computes aR and user->maxD
    user->eps = tmpeps; // restore it

    // based on maxD we can take a max-principle stable explicit step
    //dxinvsum = 1.0 / ((1.0 / user->dx) + (1.0 / user->dy));  <- replaces dx/2
    ierr = MPI_Allreduce(&user->maxD,&maxD,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);
    deltat = (user->dx / 2.0) / (2.0 * maxD);

    // take step aH, which is in-place since we have residual
    // mu = deltat / (user->dx * user->dy);
    mu = deltat / (user->dx * user->dx);
    for (k=info.ys; k<info.ys+info.ym; k++) {
        for (j=info.xs; j<info.xs+info.xm; j++) {
            aH[k][j] = PetscMax(0.0, aH[k][j] - mu * aR[k][j]);
        }
    }

    // clean up
    ierr = DMDAVecRestoreArray(user->da, H, &aH);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, R, &aR);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, Hloc, &aHloc);CHKERRQ(ierr);
    ierr = VecDestroy(&Hloc); CHKERRQ(ierr);
    ierr = VecDestroy(&R); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


