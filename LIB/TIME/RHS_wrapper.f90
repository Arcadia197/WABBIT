!> \brief wrapper for RHS call in time step function, computes RHS in work array
!! (inplace)
!!
!! calls RHS depending on physics
!!
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!**********************************************************************************************

subroutine RHS_wrapper(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID, stage)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_rhs  (:, :, :, :, :)       !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)      !> hvy_mask: the mask function, usx,usy,usz,color,sponge
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik), intent(in), optional        :: stage

    real(kind=rk), dimension(3)         :: volume_int                   !> global integral
    real(kind=rk), dimension(3)         :: dx, x0                       !> spacing and origin of a block
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, hvy_id  ! loop variables
    integer(kind=ik)                    :: g, hvy_id_mask, Bs(1:3), grhs, Neqn_RHS, Jmin_active, Jmax_active, stage_set
    integer(kind=2)                     :: n_domain(1:3)
    real(kind=rk)                       :: t0, t1, norm_b(1:1)
    real(kind=rk), save                 :: b_norm_stage1 = 0.0_rk   !< max|D S F| of the classical stage-1 solve of this step
    logical                             :: is_equidistant, incremental, start_zero
    real(kind=rk), allocatable          :: norm_ref                 !< reference norm of the relative stop test (allocated: increment solves)
    integer(kind=ik)                    :: i_p1                     !< slot of hvy_tmp that keeps p_1 for stage 4, 0: none
    character(len=cshort)               :: order_predictor_saved

    Neqn_RHS = params%n_eqn_rhs
    Bs   = params%Bs
    g    = params%g
    grhs = params%g_rhs
    t0       = MPI_wtime()
    n_domain = 0
    x0       = 0.0_rk
    dx       = 0.0_rk
    hvy_ID   = 1

    Jmin_active = minActiveLevel_tree(tree_ID_flow)
    Jmax_active = maxActiveLevel_tree(tree_ID_flow)
    is_equidistant = (Jmin_active == Jmax_active)

    ! assume this is stage 1 if we don't know. This is just important for optimizations for the Poisson solver
    stage_set = 1
    if (present(stage)) stage_set = stage
 
    !-------------------------------------------------------------------------
    ! CoarseExtension update of input data
    !-------------------------------------------------------------------------
    ! we assume data are sync'ed on call
    ! call coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
    ! hvy_n(tree_ID), inputDataSynced=.true. )

    !-------------------------------------------------------------------------
    ! create mask function at current time
    !-------------------------------------------------------------------------
    t1 = MPI_wtime()
    call createMask_tree(params, time, hvy_mask, hvy_tmp)
    call toc( "RHS_wrapper::createMask_tree", 31, MPI_wtime()-t1 )


    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    t1 = MPI_wtime()
    ! performs initializations in the RHS module, such as resetting integrals
    ! for this stage, just pass any block (result does not depend on block), hvy_id=1 and set x0=dx=0
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,1), g, x0, dx, &
         hvy_rhs(:,:,:,:,1), hvy_mask(:,:,:,:,1), "init_stage" )

    !-------------------------------------------------------------------------
    ! 2nd stage: integral_stage. (called for all blocks)
    !-------------------------------------------------------------------------
    ! For some RHS, the eqn depend not only on local, block based qtys, such as
    ! the state vector, but also on the entire grid, for example to compute a
    ! global forcing term (e.g. in FSI the forces on bodies). As the physics
    ! modules cannot see the grid, (they only see blocks), in order to encapsulate
    ! them nicer, two RHS stages have to be defined: integral / local stage.
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
        ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
        ! Note: hvy_mask is not used if in this case by the RHS routines...
        hvy_id_mask = hvy_id
        if (size(hvy_mask,5) == 1) hvy_id_mask = 1

        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, x0, dx,&
        hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id_mask), "integral_stage", n_domain )
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used ror example for MPI_REDUCES
    ! for this stage, just pass any block (result does not depend on block), hvy_id=1 and set x0=dx=0
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, 1), g, x0, dx, &
    hvy_rhs(:,:,:,:,1), hvy_mask(:,:,:,:,1), "post_stage" )
    call toc( "RHS_wrapper::integral-stage", 32, MPI_wtime()-t1 )


    !-------------------------------------------------------------------------
    ! 4th stage: local evaluation of RHS on all blocks (called for all blocks)
    !-------------------------------------------------------------------------
    ! the second stage then is what you would usually do: evaluate local differential
    ! operators etc.

    t1 = MPI_wtime()
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
        ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
        ! Note: hvy_mask is not used if in this case by the RHS routines...
        hvy_id_mask = hvy_id
        if (size(hvy_mask,5) == 1) hvy_id_mask = 1

        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, &
        x0, dx, hvy_rhs(:,:,:,:, hvy_id), hvy_mask(:,:,:,:,hvy_id_mask), "local_stage", n_domain)
    enddo
    call toc( "RHS_wrapper::local-stage", 33, MPI_wtime()-t1 )

    ! For purely hyperbolic equations, at this stage the RHS is completely ready.
    ! However, for NSPP we have to actually solve the pressure-poisson equation to get the correct pressure gradient term.
    ! This has to be done at exactly this position. However, we will need to sync some variables, which you usually would love to avoid in the RHS_wrapper
    if (params%PDE_type == "mixed_hyperbolic_parabolic") then

        ! incremental pressure (poisson_incremental): from stage 2 on, subtract the gradient of the pressure p_old of the previous
        ! stage from the RHS, so that the Poisson problem below gives only its increment dp. p_old is in hvy_block, as the
        ! time stepper leaves the pressure alone. Stage 1 solves classically, as the grid changes between the time steps
        incremental = params%poisson_incremental .and. stage_set > 1
        i_p1 = stage1_pressure_slot(params, size(hvy_tmp,4))
        if (incremental) then
            call stage4_start_pressure(params, stage_set, i_p1, hvy_block, hvy_tmp, tree_ID)
            ! the time stepper syncs only 1:n_eqn_rhs. Sync with the Poisson predictor, as the solver does for the pressure it
            ! returns, so that G p_old + G dp is the gradient of p_old + dp also in the interface strips; linear in
            ! theory, but without this sync the fields change at rounding level (1e-16, adaptive 2D), so it is kept
            order_predictor_saved = params%order_predictor
            params%order_predictor = params%poisson_order_predictor
            call sync_ghosts_tree( params, hvy_block(:,:,:,params%dim+1:params%dim+1,:), tree_ID )
            params%order_predictor = order_predictor_saved
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
                call get_block_spacing_origin( params, lgt_id, x0, dx )
                if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
                hvy_id_mask = hvy_id
                if (size(hvy_mask,5) == 1) hvy_id_mask = 1
                call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
                    hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id_mask), "pressure_gradient_stage", n_domain )
            enddo
        endif

        ! sync RHS
        t1 = MPI_wtime()
        call sync_ghosts_RHS_tree( params, hvy_rhs(:,:,:,1:Neqn_RHS, :), tree_ID, g_minus=grhs, g_plus=grhs )
        call toc( "RHS_wrapper::sync-RHS", 10035, MPI_wtime()-t1 )

        ! ! debug save
        ! call saveHDF5_tree("RHSx_0.h5", 0.0_rk, 0, 1, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)
        ! call saveHDF5_tree("RHSy_0.h5", 0.0_rk, 0, 2, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)

        ! divergence of RHS
        t1 = MPI_wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            if ( .not. All(params%periodic_BC) ) then
                ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
            endif
            ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
            ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
            ! Note: hvy_mask is not used if in this case by the RHS routines...
            hvy_id_mask = hvy_id
            if (size(hvy_mask,5) == 1) hvy_id_mask = 1
            ! RHS_meta will compute the divergence and store it in the first component of hvy_block

            call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
             hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id_mask), "divergence_stage", n_domain )
        enddo
        call toc( "RHS_wrapper::divergence-stage", 10036, MPI_wtime()-t1 )

        ! ! debug save
        ! call saveHDF5_tree("divRHS_0.h5", 0.0_rk, 0, 1, params, hvy_block(:,:,:,:,:), tree_ID, no_sync=.false.)

        ! incremental pressure: the increments of stages 2+ get the stop target of the classical stage-1 solve,
        ! tol_rel * max|D S F|, i.e. their relative test refers to the norm of the stage-1 RHS (their own RHS is much smaller,
        ! so tol_rel relative to it would be a far stricter target than the classical solve has)
        if (params%poisson_incremental .and. stage_set == 1) then
            call componentWiseNorm_tree(params, hvy_block(:,:,:,1:1,:), tree_ID, "Linfty", norm_b, threshold_state_vector=.false.)
            b_norm_stage1 = norm_b(1)
        endif
        if (incremental) norm_ref = b_norm_stage1

        ! MultiGrid solve, we sync before and after in multigrid_solve
        ! Logic here: We can reuse the latest solution as initialization if the grid is equidistant and we are sure it is the same as before
        ! (the increment dp always starts from zero)
        start_zero = incremental .or. .not.(is_equidistant .and. (stage_set > 1 .or. (stage_set == 1 .and. .not. params%adapt_tree)))
        t1 = MPI_wtime()
        if (params%poisson_exact_stages) then
            ! exact projection of the RHS: solve D G p = D rhs with the stage's own D (RHS sync, NSPP stencils), so that
            ! the pressure gradient removes the divergence of the RHS up to the solver residual. All vectors of that
            ! solve must lie in one array that the V-cycle may load-balance: hvy_tmp, which is free during the stage.
            ! b and p are copied there, and p back
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                hvy_tmp(:,:,:,1,hvy_id) = hvy_block(:,:,:,1,hvy_id)
                hvy_tmp(:,:,:,2,hvy_id) = hvy_block(:,:,:,params%dim+1,hvy_id)
            enddo
            ! incremental: dp starts from zero, the old p stays in hvy_block and is temporarily stored in the then free b slot. The
            ! last slot of hvy_tmp may keep p_1 (i_p1), then it is no Krylov slot. Every stage solve does at least one V-cycle
            call projection_solve(params, hvy_tmp, 2, 1, 3, (/ (k, k = 4, size(hvy_tmp,4) - merge(1, 0, i_p1 > 0)) /), &
                tree_ID, params_nspp%discretization, sync_rhs=.true., init_0=start_zero, time=time, &
                solve_id=stage_set, n_vc_min=1, norm_ref=norm_ref)
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                if (incremental) hvy_tmp(:,:,:,1,hvy_id) = hvy_block(:,:,:,params%dim+1,hvy_id)
                hvy_block(:,:,:,params%dim+1,hvy_id) = hvy_tmp(:,:,:,2,hvy_id)
            enddo
        else
            ! incremental: dp is solved in the pressure slot, so the old p is temporarily stored in hvy_tmp (slot 1 is free during the stage;
            ! multigrid_solve may load-balance hvy_block but restores its block order, so the temporarily stored p stays aligned)
            if (incremental) then
                do k = 1, hvy_n(tree_ID)
                    hvy_id = hvy_active(k, tree_ID)
                    hvy_tmp(:,:,:,1,hvy_id) = hvy_block(:,:,:,params%dim+1,hvy_id)
                enddo
            endif
            call multigrid_solve(params, hvy_block(:,:,:,params%dim+1:params%dim+1,:), hvy_block(:,:,:,1:1,:), hvy_block(:,:,:,2:2,:), tree_ID, init_0=start_zero, verbose=.false., hvy_full=hvy_block, time=time, solve_id=stage_set, norm_ref=norm_ref)
        endif

        ! ! debug save
        ! call saveHDF5_tree("phi_0.h5", 0.0_rk, 0, params%dim+1, params, hvy_block(:,:,:,:,:), tree_ID, no_sync=.false.)

        ! update RHS with pressure gradient
        t1 = MPI_wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            if ( .not. All(params%periodic_BC) ) then
                ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
            endif
            ! the hvy_mask array is allocated even if the mask is not used, it has then the size (1,1,1,1,1)
            ! (a single point). Therefore, pay attention not to pass hvy_mask(:,:,:,:,hvy_id) with hvy_id>1.
            ! Note: hvy_mask is not used if in this case by the RHS routines...
            hvy_id_mask = hvy_id
            if (size(hvy_mask,5) == 1) hvy_id_mask = 1

            ! RHS_meta will assume the pressure to be in component dim+1, it then will update the hvy_rhs with the pressure gradient term
            call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
             hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id_mask), "pressure_gradient_stage", n_domain )

            ! incremental pressure: the RHS now holds F - G p_old - G dp, the pressure becomes p = p_old + dp (ghost nodes
            ! included, both are synced) and is the old pressure of the next stage
            if (incremental) hvy_block(:,:,:,params%dim+1,hvy_id) = hvy_block(:,:,:,params%dim+1,hvy_id) + hvy_tmp(:,:,:,1,hvy_id)
        enddo
        call stage1_store_pressure(params, stage_set, i_p1, hvy_block, hvy_tmp, tree_ID)
        call toc( "RHS_wrapper::NSPP-pressure-poisson-solver", 10037, MPI_wtime()-t1 )

        ! ! debug save
        ! call saveHDF5_tree("RHS2x_0.h5", 0.0_rk, 0, 1, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)
        ! call saveHDF5_tree("RHS2y_0.h5", 0.0_rk, 0, 2, params, hvy_rhs(:,:,:,:,:), tree_ID, no_sync=.false.)

    endif

    ! This subroutine is actually a hack. The RHS of the ACM eqn can be equipped with a pressure gradient term that is supposed 
    ! to keep the mean flow constant - but that seems to work only approximately. I do not know where the 
    ! difference comes from. As the mean pressure forcing does not perfectly set the mean flow constant,
    ! we remove in this function the part that the forcing did not manage to remove completely. As it was now useless, we removed
    ! the mean pressure gradient from the RHS_ACM to simplify the code.
    ! Because it is a tree-level routine, it cannot be merged easily in the ACM module.
    ! -TE
    if (params_acm%use_channel_forcing) then
        call ACM_remove_channel_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
    endif
    if (params_acm%HIT_linear_forcing) then
        call HIT_remove_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
    endif

    call toc( "RHS_wrapper (TOTAL)", 30, MPI_wtime()-t0 )
end subroutine RHS_wrapper



!!! The following is a very special optimization - for incremental pressure,
!   we compute dp = p_new - p_old. It acts like a warm-start, and this is more exact if we use the pressure of the previous stage at the time of !   the current stage, which sometimes can be done via extrapolation. For classical RK4, the stage times are 0, 1/2, 1/2, 1 and we can extrapolate the fourth stage pressure: p_old = 2 p_3 - p_1

!> Slot of hvy_tmp that keeps the stage-1 pressure p_1 for the start pressure of stage 4 (see stage4_start_pressure), 0 if
!! there is none (then the optimization is silently off). hvy_tmp is free between the stages (only the mask generation writes it,
!! and only with penalization and force_maxlevel_dealiasing, then no slot); p_1 takes its last slot if the stage solve does not need it: the
!! exact projection uses the first 2m of the free slots 4:nslots, m = max(1, min(poisson_exact_m, n_free/2)), so the last slot
!! is free if removing it leaves m unchanged (the L solver uses none of them). Only for the stage times 0, 1/2, 1/2, 1 of
!! the classical RK4, which the time extrapolation of stage 4 assumes
integer(kind=ik) function stage1_pressure_slot(params, nslots)
    implicit none
    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(in)    :: nslots
    integer(kind=ik)                :: n_free

    stage1_pressure_slot = 0
    if (.not. params%poisson_incremental .or. (params%penalization .and. params%force_maxlevel_dealiasing)) return
    if (size(params%butcher_tableau,1) /= 5) return
    if (any(abs(params%butcher_tableau(1:4,1) - (/ 0.0_rk, 0.5_rk, 0.5_rk, 1.0_rk /)) > 1.0e-12_rk)) return
    n_free = nslots - 3
    if (params%poisson_exact_stages) then
        if (n_free < 3) return
        if (min(params%poisson_exact_m, (n_free-1)/2) < min(params%poisson_exact_m, n_free/2)) return
    endif
    stage1_pressure_slot = nslots
end function stage1_pressure_slot


!> Incremental pressure, stage 4: p_old = 2 p_3 - p_1, the linear extrapolation to the stage time 1 from the stage times 1/2
!! (p_3, the pressure of the previous stage, in hvy_block) and 0 (p_1, kept in slot i_p1 of hvy_tmp), instead of p_3 alone.
!! In every other stage, and without slot (i_p1 = 0), nothing is done: p_old stays the pressure of the previous stage
subroutine stage4_start_pressure(params, stage, i_p1, hvy_block, hvy_tmp, tree_ID)
    implicit none
    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(in)    :: stage, i_p1, tree_ID
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(in)       :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                :: k, hvy_id

    if (stage /= 4 .or. i_p1 == 0) return
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        hvy_block(:,:,:,params%dim+1,hvy_id) = 2.0_rk * hvy_block(:,:,:,params%dim+1,hvy_id) - hvy_tmp(:,:,:,i_p1,hvy_id)
    enddo
end subroutine stage4_start_pressure


!> Incremental pressure, stage 1: keep p_1 (ghost nodes included) in slot i_p1 of hvy_tmp for stage4_start_pressure.
!! In every other stage, and without slot (i_p1 = 0), nothing is done
subroutine stage1_store_pressure(params, stage, i_p1, hvy_block, hvy_tmp, tree_ID)
    implicit none
    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(in)    :: stage, i_p1, tree_ID
    real(kind=rk), intent(in)       :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout)    :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                :: k, hvy_id

    if (stage /= 1 .or. i_p1 == 0) return
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        hvy_tmp(:,:,:,i_p1,hvy_id) = hvy_block(:,:,:,params%dim+1,hvy_id)
    enddo
end subroutine stage1_store_pressure



subroutine ACM_remove_channel_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_rhs  (:, :, :, :, :)     !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)      !> hvy_mask: the mask function, usx,usy,usz,color,sponge
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk)                :: rhs_mean, dV, V_channel, y, safety_wall, u_bulk, tol, weight, y_lower, y_upper
    integer(kind=ik)             :: ix, iy, iz, k, hvy_id, mpierr, lgt_id, g, Bs(1:3)
    real(kind=rk), dimension(3)  :: dx, x0

    ! This subroutine is actually a hack. The RHS of the ACM eqn can be equipped with a pressure gradient term that is supposed 
    ! to keep the mean flow constant - but that seems to work only approximately. I do not know where the 
    ! difference comes from. As the mean pressure forcing does not perfectly set the mean flow constant,
    ! we remove in this function the part that the forcing did not manage to remove completely. As it was now useless, we removed
    ! the mean pressure gradient from the RHS_ACM to simplify the code.
    ! Because it is a tree-level routine, it cannot be merged easily in the ACM module.
    ! -TE
    g = params%g 
    Bs = params%Bs

    ! --- Step1: compute the mean rhs
    u_bulk = 0.0_rk
    tol = 1.0e-12_rk
    y_lower = params_acm%h_channel
    y_upper = params_acm%domain_size(2) - params_acm%h_channel
    
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )
        ! volume element
        dV = product(dx)

        ! safety distance to avoid round-off problems at the wall location
        safety_wall = dx(2) * 0.1_rk

        do iy = g+1, Bs(2)+g
            y = x0(2) + dble(iy-(g+1)) * dx(2)

            if ((y >= y_lower - tol) .and. (y <= y_upper + tol)) then

                        ! here we use the trapez-rule to solve the ingetral for u_bulk
                        ! upper and lower channel wall: 1/2
                        ! inner intervals: 1
                        if (abs(y - y_lower) <= tol .or. abs(y - y_upper) <= tol) then
                            weight = 0.5_rk
                        else
                            weight = 1.0_rk
                        endif

                        u_bulk = u_bulk + weight * sum(hvy_rhs(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 1, hvy_id)) * dV
            endif
        enddo


     enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE, u_bulk, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
    
    ! volume of fluid in the channel with walls
    V_channel = params_acm%domain_size(1)*params_acm%domain_size(3)*(params_acm%domain_size(2)-2.0_rk*params_acm%h_channel)
    u_bulk = u_bulk/V_channel

    ! --- Step 2: remove it from RHS, so that the zero mode is not drifting
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )
            
        do iy = g+1, Bs(2)+g
                y = x0(2) + dble(iy-(g+1)) * dx(2)

                if ((y >= y_lower - tol) .and. (y <= y_upper + tol)) then
                    hvy_rhs(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 1, hvy_id) = hvy_rhs(g+1:Bs(1)+g, iy, g+1:Bs(3)+g, 1, hvy_id) - u_bulk
                endif
        enddo

    enddo

end subroutine


subroutine HIT_remove_meanrhs(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, tree_ID)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_rhs  (:, :, :, :, :)     !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)      !> hvy_mask: the mask function, usx,usy,usz,color,sponge
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    real(kind=rk)                :: rhs_mean, dV, V_domain, u_RMS
    integer(kind=ik)             :: ix, iy, iz, k, hvy_id, mpierr, lgt_id, g, Bs(1:3)
    real(kind=rk), dimension(3)  :: dx, x0

    ! This subroutine is actually a hack. It enforces zero mean flow for HIT simulations with linear forcing.
    ! Because it is a tree-level routine, it cannot be merged easily in the ACM module.
    ! -TE
    g = params%g 
    Bs = params%Bs

    ! --- mean flow was already computed in RHS_ACM, we just need to remove it here
    
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! remove mean flow
        hvy_block(:,:,:,1,hvy_id) = hvy_block(:,:,:,1,hvy_id) - params_acm%mean_flow(1)
        hvy_block(:,:,:,2,hvy_id) = hvy_block(:,:,:,2,hvy_id) - params_acm%mean_flow(2)
        if (params%dim == 3) hvy_block(:,:,:,3,hvy_id) = hvy_block(:,:,:,3,hvy_id) - params_acm%mean_flow(3)
    enddo

end subroutine HIT_remove_meanrhs