subroutine threshold_block( params, u, thresholding_component, refinement_status, norm, level, input_is_WD, indices, eps, verbose_check)
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array
    !! When input_is_WD is set they are expected to be already wavelet decomposed in Spaghetti-ordering
    real(kind=rk), intent(inout)        :: u(:, :, :, :)
    !> it can be useful not to consider all components for thresholding here.
    !! e.g. to work only on the pressure or vorticity.
    logical, intent(in)                 :: thresholding_component(:)
    !> main output of this routine is the new satus
    integer(kind=ik), intent(out)       :: refinement_status
    !> If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
    !! we pass the level to this routine
    integer(kind=ik), intent(in)        :: level
    logical, intent(in)                 :: input_is_WD                       !< flag if hvy_block is already wavelet decomposed
    real(kind=rk), intent(inout)        :: norm( size(u,4) )
    !> if different from the default eps (params%eps), you can pass a different value here. This is optional
    !! and used for example when thresholding the mask function.
    real(kind=rk), intent(in), optional :: eps
    !> Indices of patch if not the whole interior block should be tresholded, used for securityZone
    integer(kind=ik), intent(in), optional :: indices(1:2, 1:3)
    logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

    integer(kind=ik)                    :: dF, i, j, l, p, idx(2,3)
    real(kind=rk)                       :: detail( size(u,4) )
    integer(kind=ik)                    :: g, i_dim, dim, Jmax, nc
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk)                       :: eps_use
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: u_wc

    nc     = size(u, 4)
    Bs     = params%Bs
    g      = params%g
    dim    = params%dim
    Jmax   = params%Jmax
    detail = -1.0_rk

    if (allocated(u_wc)) then
        if (size(u_wc, 4) < nc) deallocate(u_wc)
    endif
    if (.not. allocated(u_wc)) allocate(u_wc(1:size(u, 1), 1:size(u, 2), 1:size(u, 3), 1:nc ) )

    ! set the indices we want to treshold
    idx(:, :) = 1
    if (present(indices)) then
        idx(:, :) = indices(:, :)
    else  ! full interior block
        idx(1, 1) = g+1
        idx(2, 1) = Bs(1)+g
        idx(1, 2) = g+1
        idx(2, 2) = Bs(2)+g
        if (dim == 3) then
            idx(1, 3) = g+1
            idx(2, 3) = Bs(3)+g
        endif
    endif

    ! we need to know if the first point is a SC or WC for the patch we check and skip it if it is a WC
    !     1 2 3 4 5 6 7 8 9 A B C
    !     G G G S W S W S W G G G
    !                 I I I
    ! 1-C - index numbering in hex format, G - ghost point, S - SC, W - WC, I - point of patch to be checked
    ! Patch I is checked, but we need to know that index 7 has a WC and should be skipped
    ! this is for parity with inflatedMallat version where SC and WC are situated on the SC indices of the spaghetti format
    ! for g=odd, the SC are on even numbers; for g=even, the SC are on odd numbers
    idx(1, 1:params%dim) = idx(1, 1:params%dim) + modulo(g + idx(1, 1:params%dim) + 1, 2)
    ! also, when the last point is a SC, the last point is only partially included but its WC have to be considered
    ! this gives problem if the last point is a SC so we need to handle this special case
    do i_dim = 1, params%dim
        if (idx(2, i_dim) /= size(u, i_dim)) then
            idx(2, i_dim) = idx(2, i_dim) + modulo(idx(2, i_dim) - idx(1, i_dim) + 1, 2)
        endif
    enddo


#ifdef DEV
    if (.not. allocated(params%GD)) call abort(1213149, "The cat is angry: Wavelet-setup not yet called?")
    if (modulo(Bs(1),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
    if (modulo(Bs(2),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
#endif

    if (.not. input_is_WD) then
        u_wc(:, :, :, 1:nc) = u(:, :, :, 1:nc)  ! Wavelet decompose full block
        call waveletDecomposition_block(params, u_wc(:, :, :, 1:nc)) ! data on u (WC/SC) now in Spaghetti order
    else
        ! copy only part we need
        u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = &
           u(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc)
    endif

    ! set sc to zero to more easily compute the maxval, use offset if first index is not a SC
    u_wc(idx(1,1):idx(2,1):2, idx(1,2):idx(2,2):2, idx(1,3):idx(2,3):2, 1:nc) = 0.0_rk

    ! for L2 norm, the cross components need to be renormalized
    ! JB ToDo: Maybe this can be done more clever, maybe not multiply but simply have several maxima values?
    if (params%eps_norm == "L2") then
        ! components W_X, W_Y, W_Z    have factor 2**(dim/2-1)
        !            W_xy, W_xz, W_yz have factor 2**(dim/2-2)
        !            W_xyz             has factor 2**(dim/2-3)
        ! We treat this by multiplying all values by 2**(dim/2) and then dividing by 2 in each direction for the WC values
        u_wc(:, :, :, 1:nc) = u_wc(:, :, :, 1:nc) * 2.0_rk**(dble(params%dim)/2.0_rk)
        u_wc(idx(1,1)+1:idx(2,1)+1:2, :, :, 1:nc) = u_wc(idx(1,1)+1:idx(2,1)+1:2, :, :, 1:nc) / 2.0_rk
        u_wc(:, idx(1,2)+1:idx(2,2)+1:2, :, 1:nc) = u_wc(:, idx(1,2)+1:idx(2,2)+1:2, :, 1:nc) / 2.0_rk
        if (params%dim == 3) then
            u_wc(:, :, idx(1,3)+1:idx(2,3)+1:2, 1:nc) = u_wc(:, :, idx(1,3)+1:idx(2,3)+1:2, 1:nc) / 2.0_rk
        endif
    endif

    do p = 1, nc
        ! if all details are smaller than C_eps, we can coarsen, check interior WC only
        detail(p) = maxval( abs(u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p)) )
    enddo

    detail(1:nc) = detail(1:nc) / norm(1:nc)

    ! We could disable detail checking for qtys we do not want to consider,
    ! but this is more work and selective thresholding is rarely used
    do p = 1, nc
        if (.not. thresholding_component(p)) detail(p) = 0.0_rk
    enddo

    ! also wir brauchen einen scale(level)- dependent threshold, d.h. \epsilon_j
    ! zudem ist dieser abhaengig von der raum dimension d.
    !
    ! Fuer die L^2 normalisierung (mit wavelets welche in der L^\infty norm normalisiert sind) haben wir
    !
    ! \epsilon_j = 2^{-jd/2} \epsilon
    !
    ! d.h. der threshold wird kleiner auf kleinen skalen.
    !
    ! Fuer die vorticity (anstatt der velocity) kommt nochmal ein faktor 2^{-j} dazu, d.h.
    !
    ! \epsilon_j = 2^{-j(d+2)/2} \epsilon
    !
    ! Zum testen waere es gut in 1d oder 2d zu pruefen, ob die L^2 norm von u - u_\epsilon
    ! linear mit epsilon abnimmt, das gleiche koennte man auch fuer H^1 (philipp koennte dies doch mal ausprobieren?).
    !
    ! fuer CVS brauchen wir dann noch \epsilon was von Z (der enstrophy) und der feinsten
    ! aufloesung abhaengt. fuer L^2 normalisierte wavelets ist
    ! der threshold:
    !
    ! \epsilon = \sqrt{2/3 \sigma^2 \ln N}
    !
    ! wobei \sigma^2 die varianz (= 2 Z) der incoh. vorticity ist.
    ! typischerweise erhaelt man diese mit 1-3 iterationen.
    ! als ersten schritt koennen wir einfach Z der totalen stroemung nehmen.
    ! N ist die maximale aufloesung, typicherweise 2^{d J}.
    !

    ! default thresholding level is the one in the parameter struct
    eps_use = params%eps
    ! but if we pass another one, use that.
    if (present(eps)) eps_use = eps

    ! write(*, '("Detail ", es8.1, " eps ", es8.1)') detail(1), eps_use

    select case(params%eps_norm)
    case ("Linfty")
        ! do nothing, our wavelets are normalized in L_infty norm by default, hence
        ! a simple threshold controls this norm
        eps_use = eps_use

    case ("L2")
        ! If we want to control the L2 norm (with wavelets that are normalized in Linfty norm)
        ! threshold has to be level dependent
        eps_use = eps_use * ( 2.0_rk**(+dble((level-Jmax)*params%dim)/2.0_rk) )

    case ("H1")
        ! H1 norm mimicks filtering of vorticity
        eps_use = eps_use * ( 2**(-level*(params%dim+2.0_rk)*0.5_rk) )

    case default
        call abort(20022811, "ERROR:threshold_block.f90:Unknown wavelet normalization!")

    end select

    ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
    ! the block is tagged as "wants to coarsen" by setting the tag -1
    ! note gradedness and completeness may prevent it from actually going through with that
    if ( maxval(detail) < eps_use) then
        ! coarsen block, -1
        refinement_status = -1
    else
        refinement_status = 0
    end if

    if (present(verbose_check) .and. any(detail(:) > eps_use / 2.0_rk)) then
        write(*, '(A, es10.3, A, i2, A, 10(es10.3, 2x))') "Eps: ", eps_use, ", Ref stat: ", refinement_status, ", Details: ", detail(1:nc)
    endif
end subroutine threshold_block
