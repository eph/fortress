module fortress_smc_particles_t

  use, intrinsic :: iso_fortran_env, only: wp => real64

  use json_module, only: json_core, json_value
  use fortress_particles_t, only: fortress_particles

  implicit none

  type, extends(fortress_particles) :: fortress_smc_particles

     double precision, allocatable :: loglh(:), loglhold(:), prior(:)

   contains
     procedure :: write_json

  end type fortress_smc_particles

  interface fortress_smc_particles
     module procedure new_fortress_smc_particles
  end interface fortress_smc_particles

contains     
  type(fortress_smc_particles) function new_fortress_smc_particles(npart, nvars) result(p)

    integer, optional, intent(in) :: npart, nvars

    if (present(npart)) p%npart = npart
    if (present(nvars)) p%nvars = nvars

    allocate(p%particles(p%nvars, p%npart), p%weights(p%npart), &
         p%loglh(p%npart), p%loglhold(p%npart), p%prior(p%npart))

    p%loglh = 0.0_wp
    p%loglhold = 0.0_wp
    p%prior = 0.0_wp
    p%weights = 1.0_wp / p%npart
    p%particles = 0.0_wp
    
  end function new_fortress_smc_particles
    
  subroutine write_json(self, json_node)

    class(fortress_smc_particles) :: self
    type(json_value), pointer, intent(inout) :: json_node
    type(json_core) :: json
    character(len=3) :: varname

    integer :: j

    call json%add(json_node, 'weights', self%weights)
    call json%add(json_node, 'loglh', self%loglh)
    call json%add(json_node, 'prior', self%prior)

    do j = 1, self%nvars
       write(varname, '(I3.3)') j
       call json%add(json_node, 'var.'//trim(varname), self%particles(j,:))
    end do

  end subroutine write_json


  
  subroutine cleanup(p) 
    type(fortress_smc_particles) :: p

    deallocate(p%particles, p%weights, p%loglh, p%loglhold, p%prior)

  end subroutine cleanup


end module fortress_smc_particles_t
