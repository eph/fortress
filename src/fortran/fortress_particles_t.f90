module fortress_particles_t

  use, intrinsic :: iso_fortran_env, only: wp => real64
  use json_module, only: json_core, json_value


  implicit none

  type fortress_particles

     integer :: npart = 1000
     integer :: nvars = 1

     double precision, allocatable :: particles(:,:), weights(:)

   contains


     procedure :: mean
     procedure :: mean_and_variance
     procedure :: normalize_weights
     procedure :: ESS
     procedure :: describe 
     procedure :: systematic_resampling
     procedure :: free
     procedure :: write
     procedure :: write_json
     !final :: cleanup 

  end type fortress_particles


  interface fortress_particles
     module procedure new_fortress_particles
  end interface fortress_particles

contains

  type(fortress_particles) function new_fortress_particles(npart, nvars) result(p)

    integer, optional, intent(in) :: npart, nvars

    if (present(npart)) p%npart = npart
    if (present(nvars)) p%nvars = nvars

    allocate(p%particles(p%nvars, p%npart), p%weights(p%npart))

    p%particles = 0.0d0
    p%weights = 1.0d0 / p%npart

  end function new_fortress_particles


  subroutine normalize_weights(p, z)
    class(fortress_particles) :: p
    double precision, intent(out) :: z 

    z = sum(p%weights)

    p%weights = p%weights / z 

  end subroutine normalize_weights


  function ESS(p) result(effective_sample_size)
    class(fortress_particles) :: p

    double precision :: effective_sample_size 
    double precision :: m1, logw(p%npart)

    logw = 2.0_wp*log(p%weights)
    m1 = maxval(logw)
    
    effective_sample_size = 1.0d0 / (sum(exp(logw-m1)) * exp(m1))

  end function ESS

  subroutine describe(p) 
    class(fortress_particles) :: p

    double precision :: mu(p%nvars), std(p%nvars)

    integer :: i

    print*, 'Describing Particle Swarm' 
    print*, '# variables = ', p%nvars
    print*, '# particles = ', p%npart
    print*, 'ESS         = ', p%ESS()

    mu = p%mean()
    do i = 1, p%nvars
       print*, 'variable', i, mu(i)
    end do

  end subroutine describe

  function mean(p) result(mu)
    class(fortress_particles) :: p
    double precision :: mu(p%nvars)

    integer :: i

    do i = 1,p%nvars
       mu(i) = dot_product(p%particles(i,:), p%weights)
    end do

  end function mean

  subroutine mean_and_variance(p, mean, variance)
    class(fortress_particles) :: p
    double precision, intent(out) :: mean(p%nvars), variance(p%nvars, p%nvars)
    integer :: j

    mean = p%mean()
    variance = 0.0_wp
    do j = 1, p%npart
       call dger(p%nvars, p%nvars, p%weights(j), p%particles(:,j), 1, &
            p%particles(:,j), 1, variance, p%nvars)
    end do
    call dger(p%nvars, p%nvars, -1.0_wp, mean, 1, mean, 1, variance, p%nvars)

  end subroutine mean_and_variance

  subroutine systematic_resampling(p, randu, resample_ind)
    class(fortress_particles) :: p

    double precision, intent(in) :: randu
    integer, intent(out), optional :: resample_ind(p%npart)
    double precision :: cdf(p%npart), uu(p%npart)
    double precision :: part_rep(p%nvars, p%npart)

    integer :: i,j

    !wtold = wt/sum(wt)
    cdf(1) = p%weights(1)!wtold(1)
    do i=2,p%npart
       cdf(i) = cdf(i-1) + p%weights(i)
    end do

    uu = ( randu -1.0d0 + real( (/ (i, i=1,p%npart) /) ,8) ) / real(p%npart,8)


    j=1
    do i=1,p%npart
       ! move along the CDF
       do while (uu(i)>cdf(j))
          j=j+1
       end do
       ! shuffling
       part_rep(:,i) = p%particles(:,j)

       if (present(resample_ind)) resample_ind(i) = j
    end do

    p%particles = part_rep
    p%weights = 1.0d0 / p%npart

  end subroutine systematic_resampling

  subroutine write(self, file_prefix)
    class(fortress_particles) :: self

    character(len=200), intent(in) :: file_prefix

    character(len=200) :: state_name, weight_name

    integer :: i

    state_name = trim(adjustl(trim(file_prefix)))//"states.txt"
    weight_name = trim(adjustl(trim(file_prefix)))//"weights.txt"

    open(345, file=state_name, action='write')
    open(346, file=weight_name, action='write')

    do i = 1, self%npart
       write(345, '(100f16.8)') self%particles(:, i)
       write(346, '(100f16.8)') self%weights(i)
    end do

    close(345)
    close(346)

  end subroutine write

  subroutine write_json(self, json_node)

    class(fortress_particles) :: self
    type(json_value), pointer, intent(inout) :: json_node
    type(json_core) :: json
    character(len=3) :: varname

    integer :: j
    call json%add(json_node, 'weights', self%weights)

    do j = 1, self%nvars
       write(varname, '(I3.3)') j
       call json%add(json_node, 'var.'//trim(varname), self%particles(j,:))
    end do

  end subroutine write_json
  
  subroutine free(p) 
    class(fortress_particles) :: p

    deallocate(p%particles, p%weights)

  end subroutine free

  subroutine cleanup(p) 
    type(fortress_particles) :: p

    deallocate(p%particles, p%weights)

  end subroutine cleanup


end module fortress_particles_t
