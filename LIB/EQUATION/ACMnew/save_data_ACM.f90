!-----------------------------------------------------------------------------
! save data. Since you might want to save derived data, such as the vorticity,
! the divergence etc., which are not in your state vector, this routine has to
! copy and compute what you want to save to the work array.
!
! In the main code, save_fields than saves the first N_fields_saved components of the
! work array to file.
!
! NOTE that as we have way more work arrays than actual state variables (typically
! for a RK4 that would be >= 4*dim), you can compute a lot of stuff, if you want to.
!-----------------------------------------------------------------------------
subroutine PREPARE_SAVE_DATA_ACM( time, u, g, x0, dx, work, mask, n_domain )
    implicit none
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! output in work array.
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(inout) :: mask(1:,1:,1:,1:)

    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The n_domain tells you if the local field is adjacent to a domain boundary:
    ! n_domain(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    ! NOTE: ACM only supports symmetry BC for the moment (which is handled by wabbit and not ACM)
    integer(kind=2), intent(in) :: n_domain(3)

    ! local variables
    integer(kind=ik)  :: neqn, nwork, k, iscalar, iaverage
    integer(kind=ik), dimension(3) :: Bs
    character(len=cshort) :: name

    if (.not. params_acm%initialized) write(*,*) "WARNING: PREPARE_SAVE_DATA_ACM called but ACM not initialized"
    ! number of state variables
    neqn = size(u,4)
    ! number of available work array slots
    nwork = size(work,4)


    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    if (params_acm%geometry == "Insect" ) then
        call Update_Insect(time, Insect)
    endif



    do k = 1, size(params_acm%names,1)
        name = params_acm%names(k)
        select case(name)
        case('ux', 'Ux', 'UX')
            ! copy state vector
            work(:,:,:,k) = u(:,:,:,1)

        case('uy', 'Uy', 'UY')
            ! copy state vector
            work(:,:,:,k) = u(:,:,:,2)

        case('uz', 'Uz', 'UZ')
            ! copy state vector
            work(:,:,:,k) = u(:,:,:,3)

        case('p', 'P')
            ! copy state vector (do not use 4 but rather neq for 2D runs, where p=3rd component)
            work(:,:,:,k) = u(:,:,:,params_acm%dim+1)

        case('vor', 'vort', 'Vor', 'Vort', 'vorticity', 'Vorticity')
            if (size(work,4) - k < 2) then
                call abort(19101810,"ACM: Not enough space to compute vorticity, put vorticity computation as first save variables or put atleast two other save variables afterwards. Works only in 2D (known bug)")
            endif
            if (params_acm%dim /= 2) then
                call abort(19101811,"ACM: storing scalar vor is not possible in 3D - use any of 'vorx' 'vory' 'vorz' 'vorabs' or compute in post (known bug)")
            endif

            ! vorticity
            call compute_vorticity(u(:,:,:,1:3), &
            dx, Bs, g, params_acm%discretization, work(:,:,:,k:k+2))

        case('vorx', 'Vorx', 'VorX', 'vory', 'Vory', 'VorY', 'vorz', 'Vorz', 'VorZ', 'vorabs', 'Vorabs', 'VorAbs', 'vor-abs', 'Vor-abs', 'Vor-Abs')
            if (size(work,4) - k < 2) then
                call abort(19101810,"ACM: Not enough space to compute vorticity, put vorticity computation as first save variables or put atleast two other save variables afterwards. Works only in 3D (known bug)")
            endif
            if (params_acm%dim /= 3) then
                call abort(19101811,"ACM: storing vector vor is not possible in 2D - use 'vor' or compute in post (known bug)")
            endif

            ! vorticity, this effectively computes it three times for all components, but I just assume we do not save often
            call compute_vorticity(u(:,:,:,1:3), &
            dx, Bs, g, params_acm%discretization, work(:,:,:,k:k+2))
            ! for different components y,z we need to copy the desired one to the first position
            if (name == 'vory' .or. name == 'Vory' .or. name == 'VorY') then
                work(:,:,:,k) = work(:,:,:,k+1)
            elseif (name == 'vorz' .or. name == 'Vorz' .or. name == 'VorZ') then
                work(:,:,:,k) = work(:,:,:,k+2)
            elseif (name=='vorabs' .or. name=='Vorabs' .or.name=='VorAbs' .or. name=='vor-abs' .or. name=='Vor-abs' .or. name=='Vor-Abs') then
                work(:,:,:,k) = sqrt(work(:,:,:,k)**2 + work(:,:,:,k+1)**2 + work(:,:,:,k+2)**2)
            endif
        
        case('helx', 'Helx', 'HelX', 'hely', 'Hely', 'HelY', 'helz', 'Helz', 'HelZ', 'helabs', 'Helabs', 'HelAbs', 'hel-abs', 'Hel-abs', 'Hel-Abs')
            if (size(work,4) - k < 2) then
                call abort(19101810,"ACM: Not enough space to compute vorticity for helicity, put vorticity computation as first save variables or put atleast two other save variables afterwards. Works only in 3D (known bug)")
            endif
            if (params_acm%dim /= 3) then
                call abort(19101811,"ACM: Helicity is always 0 in 2D - we do not need to store that!")
            endif

            ! vorticity, this effectively computes it three times for all components, but I just assume we do not save often
            call compute_vorticity(u(:,:,:,1:3), &
            dx, Bs, g, params_acm%discretization, work(:,:,:,k:k+2))
            ! compute by velocity to get local helicity
            work(:,:,:,k:k+2) = work(:,:,:,k:k+2) * u(:,:,:,1:3)
            ! for different components y,z we need to copy the desired one to the first position
            if (name == 'hely' .or. name == 'Hely' .or. name == 'HelY') then
                work(:,:,:,k) = work(:,:,:,k+1)
            elseif (name == 'helz' .or. name == 'Helz' .or. name == 'HelZ') then
                work(:,:,:,k) = work(:,:,:,k+2)
            elseif (name=='helabs' .or. name=='Helabs' .or.name=='HelAbs' .or. name=='hel-abs' .or. name=='Hel-abs' .or. name=='Hel-Abs') then
                work(:,:,:,k) = sqrt(work(:,:,:,k)**2 + work(:,:,:,k+1)**2 + work(:,:,:,k+2)**2)
            endif

        case('div', 'divu', 'divergence')
            ! div(u)
            call compute_divergence(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,k))
        
        case('diss', 'dissipation')
            ! local dissipation rate computed over velocity weighted laplacian
            call compute_dissipation(u(:,:,:,1:params_acm%dim), dx, Bs, g, params_acm%discretization, work(:,:,:,k))
        
        case('gradPx', 'gradPy', 'gradPz', 'gradpx', 'gradpy', 'gradpz', 'gradientpx', 'gradientpy', 'gradientpz', 'gradientPx', 'gradientPy', 'gradientPz', &
            'grad-Px', 'grad-Py', 'grad-Pz', 'grad-px', 'grad-py', 'grad-pz', 'gradient-px', 'gradient-py', 'gradient-pz', 'gradient-Px', 'gradient-Py', 'gradient-Pz')
            ! Gradient of pressure, I admit the naming is confusing so I just added every option I could think of
            if (INDEX(name, 'x') > 0) then
                call compute_derivative(u(:,:,:,params_acm%dim+1), dx, Bs, g, 1, 1, params_acm%discretization, work(:,:,:,k))
            elseif (INDEX(name, 'y') > 0) then
                call compute_derivative(u(:,:,:,params_acm%dim+1), dx, Bs, g, 2, 1, params_acm%discretization, work(:,:,:,k))
            elseif (INDEX(name, 'z') > 0) then
                if (params_acm%dim == 3) then
                    call compute_derivative(u(:,:,:,params_acm%dim+1), dx, Bs, g, 3, 1, params_acm%discretization, work(:,:,:,k))
                else
                    call abort(19101812,"ACM: Gradient of p in z-direction is not defined for 2D runs")
                endif
            endif  

        case('mask')
            work(:,:,:,k) = mask(:,:,:,1)

        case('usx')
            work(:,:,:,k) = mask(:,:,:,2)

        case('usy')
            work(:,:,:,k) = mask(:,:,:,3)

        case('usz')
            work(:,:,:,k) = mask(:,:,:,4)

        case('color')
            work(:,:,:,k) = mask(:,:,:,5)

        case('sponge')
            work(:,:,:,k) = mask(:,:,:,6)

        end select


        if (name(1:6) == "scalar") then
            read( name(7:7), * ) iscalar
            work(:,:,:,k) = u(:,:,:,params_acm%dim + 1 + iscalar)
        endif
        if (name(1:14) == "timestatistics") then
            read( name(15:15), * ) iaverage
            work(:,:,:,k) = u(:,:,:,params_acm%dim + 1 + params_acm%N_scalars + iaverage)
        endif
    end do

end subroutine


!-----------------------------------------------------------------------------
! when savig to disk, WABBIT would like to know how you named you variables.
! e.g. u(:,:,:,1) is called "ux"
!
! the main routine save_fields has to know how you label the stuff you want to
! store from the work array, and this routine returns those strings
!-----------------------------------------------------------------------------
subroutine FIELD_NAMES_ACM( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=cshort), intent(out) :: name

    if (.not. params_acm%initialized) write(*,*) "WARNING: FIELD_NAMES_ACM called but ACM not initialized"

    if (allocated(params_acm%names)) then
        name = params_acm%names(N)
    else
        call abort(5554,'Something ricked')
    endif

end subroutine FIELD_NAMES_ACM
