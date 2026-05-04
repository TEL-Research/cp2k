!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief MaxwellLink coupling helpers for time-dependent dynamics.
! **************************************************************************************************
MODULE maxwelllink
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_convert_offsets_to_sizes,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_distribution_type,&
                                              dbcsr_init_p,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_trace
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE force_env_types,                 ONLY: force_env_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_result_methods,               ONLY: get_results,&
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE input_constants,                 ONLY: do_em,&
                                              nve_ensemble,&
                                              use_mom_ref_zero
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              max_line_length
   USE maxwelllink_socket,              ONLY: maxwelllink_socket_connect,&
                                              maxwelllink_socket_receive_field,&
                                              maxwelllink_socket_send_source,&
                                              maxwelllink_socket_shutdown,&
                                              maxwelllink_socket_type
   USE maxwelllink_types,               ONLY: mxl_socket_control_type
   USE md_ener_types,                   ONLY: md_ener_type
   USE md_energies,                     ONLY: md_energy
   USE md_environment_types,            ONLY: get_md_env,&
                                              md_environment_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: c_light_au
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_moments,                      ONLY: qs_moment_locop
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'maxwelllink'

   TYPE maxwelllink_type
      TYPE(maxwelllink_socket_type)          :: socket
      LOGICAL                                :: active = .FALSE.
      LOGICAL                                :: have_result = .FALSE.
      LOGICAL                                :: periodic_velocity_gauge = .FALSE.
      REAL(KIND=dp), DIMENSION(3)            :: dipole_initial = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: dipole_before_last = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: dipole_last = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: field_endpoint = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: field_midpoint = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: vector_potential = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: current_last = 0.0_dp
      REAL(KIND=dp)                          :: cell_volume = 0.0_dp
      REAL(KIND=dp)                          :: energy_last = 0.0_dp
      INTEGER                                :: dipole_history = 0
      INTEGER                                :: current_history = 0
   END TYPE maxwelllink_type

   TYPE maxwelllink_ehrenfest_type
      TYPE(maxwelllink_type)                  :: link
      TYPE(dft_control_type), POINTER         :: dft_control => NULL()
   END TYPE maxwelllink_ehrenfest_type

   PUBLIC :: maxwelllink_type,&
             maxwelllink_ehrenfest_type,&
             maxwelllink_ehrenfest_init,&
             maxwelllink_ehrenfest_receive_field,&
             maxwelllink_ehrenfest_send_source,&
             maxwelllink_ehrenfest_shutdown,&
             maxwelllink_ehrenfest_validate,&
             maxwelllink_init,&
             maxwelllink_receive_field,&
             maxwelllink_send_source,&
             maxwelllink_shutdown,&
             maxwelllink_validate_rtp_control

CONTAINS

! **************************************************************************************************
!> \brief Caches MaxwellLink RTP control for Ehrenfest MD hooks.
!> \param mxl MaxwellLink Ehrenfest runtime state.
!> \param force_env Force environment.
! **************************************************************************************************
   SUBROUTINE maxwelllink_ehrenfest_attach(mxl, force_env)
      TYPE(maxwelllink_ehrenfest_type), INTENT(INOUT)    :: mxl
      TYPE(force_env_type), POINTER                      :: force_env

      IF (.NOT. ASSOCIATED(mxl%dft_control)) &
         CALL get_qs_env(force_env%qs_env, dft_control=mxl%dft_control)

   END SUBROUTINE maxwelllink_ehrenfest_attach

! **************************************************************************************************
!> \brief Validates MaxwellLink settings for Ehrenfest MD.
!> \param mxl MaxwellLink Ehrenfest runtime state.
!> \param force_env Force environment.
!> \param ensemble MD ensemble.
! **************************************************************************************************
   SUBROUTINE maxwelllink_ehrenfest_validate(mxl, force_env, ensemble)
      TYPE(maxwelllink_ehrenfest_type), INTENT(INOUT)    :: mxl
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER, INTENT(IN)                                :: ensemble

      CALL maxwelllink_ehrenfest_attach(mxl, force_env)
      IF (mxl%dft_control%rtp_control%mxl_socket%active .AND. ensemble /= nve_ensemble) &
         CPABORT("MaxwellLink Ehrenfest coupling currently supports only the NVE ensemble.")

   END SUBROUTINE maxwelllink_ehrenfest_validate

! **************************************************************************************************
!> \brief Initializes MaxwellLink coupling for Ehrenfest MD.
!> \param mxl MaxwellLink Ehrenfest runtime state.
!> \param force_env Force environment.
!> \param potential_energy Initial CP2K potential/electronic energy.
! **************************************************************************************************
   SUBROUTINE maxwelllink_ehrenfest_init(mxl, force_env, potential_energy)
      TYPE(maxwelllink_ehrenfest_type), INTENT(INOUT)    :: mxl
      TYPE(force_env_type), POINTER                      :: force_env
      REAL(KIND=dp), INTENT(IN)                          :: potential_energy

      CALL maxwelllink_ehrenfest_attach(mxl, force_env)
      CALL maxwelllink_init(mxl%link, force_env%qs_env, mxl%dft_control%rtp_control, potential_energy)

   END SUBROUTINE maxwelllink_ehrenfest_init

! **************************************************************************************************
!> \brief Receives and installs the MaxwellLink field before an Ehrenfest MD step.
!> \param mxl MaxwellLink Ehrenfest runtime state.
!> \param force_env Force environment.
!> \param stop MaxwellLink requested termination.
! **************************************************************************************************
   SUBROUTINE maxwelllink_ehrenfest_receive_field(mxl, force_env, stop)
      TYPE(maxwelllink_ehrenfest_type), INTENT(INOUT)    :: mxl
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(OUT)                               :: stop

      stop = .FALSE.
      IF (.NOT. mxl%link%active) RETURN

      CALL maxwelllink_ehrenfest_attach(mxl, force_env)
      CALL maxwelllink_receive_field(mxl%link, force_env%qs_env, mxl%dft_control%rtp_control, &
                                     force_env%qs_env%rtp%dt, stop)

   END SUBROUTINE maxwelllink_ehrenfest_receive_field

! **************************************************************************************************
!> \brief Sends the MaxwellLink source after a converged Ehrenfest MD step.
!> \param mxl MaxwellLink Ehrenfest runtime state.
!> \param md_env MD environment.
!> \param force_env Force environment.
!> \param stop MaxwellLink requested termination.
!> \param energy_current Whether md_ener was refreshed in this call.
! **************************************************************************************************
   SUBROUTINE maxwelllink_ehrenfest_send_source(mxl, md_env, force_env, stop, energy_current)
      TYPE(maxwelllink_ehrenfest_type), INTENT(INOUT)    :: mxl
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(OUT)                               :: stop
      LOGICAL, INTENT(OUT), OPTIONAL                     :: energy_current

      TYPE(md_ener_type), POINTER                        :: md_ener

      stop = .FALSE.
      IF (PRESENT(energy_current)) energy_current = .FALSE.
      IF (.NOT. mxl%link%active) RETURN

      NULLIFY (md_ener)
      CALL get_md_env(md_env, md_ener=md_ener)
      CALL md_energy(md_env, md_ener)
      IF (PRESENT(energy_current)) energy_current = .TRUE.
      CALL maxwelllink_send_source(mxl%link, force_env%qs_env, force_env%qs_env%rtp%dt, &
                                   force_env%qs_env%sim_time, md_ener%epot, stop, &
                                   kinetic_energy=md_ener%ekin)

   END SUBROUTINE maxwelllink_ehrenfest_send_source

! **************************************************************************************************
!> \brief Shuts down MaxwellLink coupling for Ehrenfest MD.
!> \param mxl MaxwellLink Ehrenfest runtime state.
!> \param force_env Force environment.
! **************************************************************************************************
   SUBROUTINE maxwelllink_ehrenfest_shutdown(mxl, force_env)
      TYPE(maxwelllink_ehrenfest_type), INTENT(INOUT)    :: mxl
      TYPE(force_env_type), POINTER                      :: force_env

      CALL maxwelllink_shutdown(mxl%link, force_env%qs_env)
      NULLIFY (mxl%dft_control)

   END SUBROUTINE maxwelllink_ehrenfest_shutdown

! **************************************************************************************************
!> \brief Validates RT-TDDFT control settings for MaxwellLink coupling.
!> \param rtp_control RTP control settings.
! **************************************************************************************************
   SUBROUTINE maxwelllink_validate_rtp_control(rtp_control)
      TYPE(rtp_control_type), INTENT(IN)                 :: rtp_control

      IF (.NOT. rtp_control%mxl_socket%active) RETURN

      IF (rtp_control%propagator /= do_em) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET requires REAL_TIME_PROPAGATION%PROPAGATOR EM.")
      IF (rtp_control%velocity_gauge .AND. rtp_control%linear_scaling) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET periodic velocity-gauge support requires MO-based RTP.")
      IF (rtp_control%velocity_gauge .AND. rtp_control%nl_gauge_transform) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET periodic velocity-gauge current source does not support VG_COM_NL.")

   END SUBROUTINE maxwelllink_validate_rtp_control

! **************************************************************************************************
!> \brief Initializes MaxwellLink coupling.
!> \param mxl MaxwellLink runtime state.
!> \param qs_env Quickstep environment.
!> \param rtp_control RTP control settings.
!> \param energy Initial CP2K potential/electronic energy.
! **************************************************************************************************
   SUBROUTINE maxwelllink_init(mxl, qs_env, rtp_control, energy)
      TYPE(maxwelllink_type), INTENT(INOUT)              :: mxl
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rtp_control_type), INTENT(INOUT)              :: rtp_control
      REAL(KIND=dp), INTENT(IN)                          :: energy

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_init'

      CHARACTER(LEN=default_path_length)                 :: socket_host
      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: periodic_cell
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      IF (.NOT. rtp_control%mxl_socket%active) THEN
         CALL timestop(handle)
         RETURN
      END IF

      NULLIFY (cell, dft_control, para_env)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, para_env=para_env)
      periodic_cell = ANY(cell%perd(1:3) /= 0)
      IF (periodic_cell .AND. .NOT. rtp_control%velocity_gauge) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET periodic RT-TDDFT requires REAL_TIME_PROPAGATION%VELOCITY_GAUGE.")
      IF ((.NOT. periodic_cell) .AND. rtp_control%velocity_gauge) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET velocity-gauge support is currently restricted to periodic cells.")
      IF (rtp_control%propagator /= do_em) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET requires REAL_TIME_PROPAGATION%PROPAGATOR EM.")
      IF (periodic_cell .AND. ASSOCIATED(dft_control%efield_fields)) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET periodic velocity-gauge mode cannot be combined with DFT%EFIELD.")

      mxl%active = .TRUE.
      mxl%have_result = .FALSE.
      mxl%periodic_velocity_gauge = periodic_cell
      mxl%energy_last = energy
      mxl%field_endpoint = 0.0_dp
      mxl%field_midpoint = 0.0_dp
      mxl%vector_potential = c_light_au*rtp_control%vec_pot
      mxl%current_last = 0.0_dp
      mxl%current_history = 0
      mxl%cell_volume = cell%deth
      socket_host = maxwelllink_socket_host(rtp_control%mxl_socket)

      IF (para_env%is_source()) THEN
         output_unit = cp_logger_get_default_io_unit()
         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(/,T2,A)") "MAXWELLLINK SOCKET"
            IF (rtp_control%mxl_socket%unix) THEN
               WRITE (output_unit, FMT="(T3,A,A)") "UNIX socket: ", TRIM(socket_host)
            ELSE
               WRITE (output_unit, FMT="(T3,A,A)") "Host: ", TRIM(socket_host)
               WRITE (output_unit, FMT="(T3,A,I0)") "Port: ", rtp_control%mxl_socket%port
            END IF
         END IF
         CALL maxwelllink_socket_connect(mxl%socket, socket_host, rtp_control%mxl_socket%port, &
                                         rtp_control%mxl_socket%unix, rtp_control%mxl_socket%verbosity)
      END IF

      IF (mxl%periodic_velocity_gauge) THEN
         mxl%dipole_initial = 0.0_dp
         mxl%dipole_last = 0.0_dp
         CALL maxwelllink_compute_periodic_current(qs_env, mxl%current_last)
         mxl%current_history = 1
      ELSE IF (rtp_control%mxl_socket%reset_dipole) THEN
         CALL maxwelllink_compute_dipole(qs_env, mxl%dipole_initial)
      ELSE
         mxl%dipole_initial = 0.0_dp
      END IF
      IF (.NOT. mxl%periodic_velocity_gauge) THEN
         CALL maxwelllink_compute_dipole(qs_env, mxl%dipole_last)
         mxl%dipole_last = mxl%dipole_last - mxl%dipole_initial
      END IF
      mxl%dipole_before_last = mxl%dipole_last
      mxl%dipole_history = 1

      CALL timestop(handle)

   END SUBROUTINE maxwelllink_init

! **************************************************************************************************
!> \brief Receives and installs the MaxwellLink field for the next dynamics step.
!> \param mxl MaxwellLink runtime state.
!> \param qs_env Quickstep environment.
!> \param rtp_control RTP control settings.
!> \param dt Dynamics timestep in atomic units.
!> \param stop MaxwellLink requested termination.
! **************************************************************************************************
   SUBROUTINE maxwelllink_receive_field(mxl, qs_env, rtp_control, dt, stop)
      TYPE(maxwelllink_type), INTENT(INOUT)              :: mxl
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rtp_control_type), INTENT(INOUT)              :: rtp_control
      REAL(KIND=dp), INTENT(IN)                          :: dt
      LOGICAL, INTENT(OUT)                               :: stop

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_receive_field'

      INTEGER                                            :: handle, received_molecule_id
      LOGICAL                                            :: received_init
      REAL(KIND=dp)                                      :: init_dt
      REAL(KIND=dp), DIMENSION(3)                        :: field
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      stop = .FALSE.
      IF (.NOT. mxl%active) THEN
         CALL timestop(handle)
         RETURN
      END IF

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      IF (para_env%is_source()) THEN
         CALL maxwelllink_socket_receive_field(mxl%socket, mxl%have_result, field, stop, &
                                               received_init, init_dt, received_molecule_id)
      ELSE
         field = 0.0_dp
         stop = .FALSE.
         received_init = .FALSE.
         init_dt = -1.0_dp
         received_molecule_id = -1
      END IF

      CALL para_env%bcast(field)
      CALL para_env%bcast(stop)
      CALL para_env%bcast(received_init)
      CALL para_env%bcast(init_dt)
      CALL para_env%bcast(received_molecule_id)

      IF (received_init) THEN
         IF (init_dt > 0.0_dp .AND. ABS(init_dt - dt) > 1.0E-10_dp*MAX(1.0_dp, ABS(dt))) &
            CALL cp_abort(__LOCATION__, "MaxwellLink INIT dt_au does not match MOTION%MD%TIMESTEP.")
         IF (rtp_control%mxl_socket%molecule_id >= 0 .AND. &
             received_molecule_id /= rtp_control%mxl_socket%molecule_id) &
            CALL cp_abort(__LOCATION__, "MaxwellLink INIT molecule id does not match MXL_SOCKET%MOLECULE_ID.")
      END IF

      IF (.NOT. stop) THEN
         mxl%field_midpoint = field
         IF (mxl%periodic_velocity_gauge) THEN
            mxl%vector_potential = mxl%vector_potential - c_light_au*dt*mxl%field_midpoint
            rtp_control%field = mxl%field_midpoint
            rtp_control%vec_pot = mxl%vector_potential/c_light_au
         ELSE
            ! MaxwellLink sends midpoint fields; EM will average this endpoint with the previous one.
            field = 2.0_dp*mxl%field_midpoint - mxl%field_endpoint
            mxl%field_endpoint = field
            rtp_control%field = field
         END IF
      END IF
      mxl%have_result = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE maxwelllink_receive_field

! **************************************************************************************************
!> \brief Sends the source-current vector for the step that just converged.
!> \param mxl MaxwellLink runtime state.
!> \param qs_env Quickstep environment.
!> \param dt Dynamics timestep in atomic units.
!> \param sim_time Current simulation time in atomic units.
!> \param energy Current CP2K potential/electronic energy, including field coupling.
!> \param stop MaxwellLink requested termination.
!> \param kinetic_energy Optional nuclear kinetic energy for Ehrenfest dynamics.
! **************************************************************************************************
   SUBROUTINE maxwelllink_send_source(mxl, qs_env, dt, sim_time, energy, stop, kinetic_energy)
      TYPE(maxwelllink_type), INTENT(INOUT)              :: mxl
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: dt, sim_time, energy
      LOGICAL, INTENT(OUT)                               :: stop
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: kinetic_energy

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_send_source'

      CHARACTER(LEN=max_line_length)                     :: extra_json
      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: energy_field_free, kinetic_field_free, &
                                                            potential_field_free, source_time
      REAL(KIND=dp), DIMENSION(3)                        :: current_density, dipole_abs, dipole_end, &
                                                            dipole_midpoint, source
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      stop = .FALSE.
      IF (.NOT. mxl%active) THEN
         CALL timestop(handle)
         RETURN
      END IF

      kinetic_field_free = 0.0_dp
      IF (PRESENT(kinetic_energy)) kinetic_field_free = kinetic_energy
      source_time = sim_time
      IF (mxl%periodic_velocity_gauge) THEN
         CALL maxwelllink_compute_periodic_current(qs_env, source)
         IF (mxl%current_history >= 1) THEN
            dipole_end(:) = mxl%dipole_last(:) + 0.5_dp*dt*(mxl%current_last(:) + source(:))
         ELSE
            dipole_end(:) = mxl%dipole_last(:) + dt*source(:)
         END IF
         dipole_abs(:) = dipole_end(:)
         dipole_midpoint(:) = 0.5_dp*(dipole_end(:) + mxl%dipole_last(:))
         potential_field_free = energy
         energy_field_free = potential_field_free + kinetic_field_free
         IF (mxl%cell_volume > 0.0_dp) THEN
            current_density(:) = source(:)/mxl%cell_volume
         ELSE
            current_density(:) = 0.0_dp
         END IF
         CALL maxwelllink_build_json_periodic(source_time, energy_field_free, kinetic_field_free, &
                                              potential_field_free, dipole_end, dipole_midpoint, source, &
                                              current_density, mxl%field_midpoint, mxl%vector_potential, &
                                              extra_json)
      ELSE
         CALL maxwelllink_compute_dipole(qs_env, dipole_abs)
         dipole_end = dipole_abs - mxl%dipole_initial
         IF (mxl%dipole_history >= 2) THEN
            source(:) = (3.0_dp*dipole_end(:) - 4.0_dp*mxl%dipole_last(:) + mxl%dipole_before_last(:))/(2.0_dp*dt)
         ELSE
            source(:) = (dipole_end(:) - mxl%dipole_last(:))/dt
         END IF
         dipole_midpoint(:) = 0.5_dp*(dipole_end(:) + mxl%dipole_last(:))
         potential_field_free = energy + DOT_PRODUCT(dipole_abs, mxl%field_endpoint)
         energy_field_free = potential_field_free + kinetic_field_free
         CALL maxwelllink_build_json(source_time, energy_field_free, kinetic_field_free, &
                                     potential_field_free, dipole_end, dipole_midpoint, extra_json)
      END IF

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      IF (para_env%is_source()) THEN
         CALL maxwelllink_socket_send_source(mxl%socket, energy_field_free, source, TRIM(extra_json), stop)
      ELSE
         stop = .FALSE.
      END IF
      CALL para_env%bcast(stop)

      mxl%have_result = .FALSE.
      mxl%dipole_before_last = mxl%dipole_last
      mxl%dipole_last = dipole_end
      mxl%dipole_history = MIN(mxl%dipole_history + 1, 2)
      IF (mxl%periodic_velocity_gauge) THEN
         mxl%current_last = source
         mxl%current_history = 1
      END IF
      mxl%energy_last = energy_field_free

      CALL timestop(handle)

   END SUBROUTINE maxwelllink_send_source

! **************************************************************************************************
!> \brief Shuts down MaxwellLink coupling.
!> \param mxl MaxwellLink runtime state.
!> \param qs_env Quickstep environment.
! **************************************************************************************************
   SUBROUTINE maxwelllink_shutdown(mxl, qs_env)
      TYPE(maxwelllink_type), INTENT(INOUT)              :: mxl
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(mp_para_env_type), POINTER                    :: para_env

      IF (.NOT. mxl%active) RETURN

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      IF (para_env%is_source()) CALL maxwelllink_socket_shutdown(mxl%socket)
      mxl%active = .FALSE.
      mxl%have_result = .FALSE.

   END SUBROUTINE maxwelllink_shutdown

! **************************************************************************************************
!> \brief Computes the non-periodic dipole moment in atomic units.
!> \param qs_env Quickstep environment.
!> \param dipole Dipole moment in atomic units.
! **************************************************************************************************
   SUBROUTINE maxwelllink_compute_dipole(qs_env, dipole)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: dipole

      CHARACTER(LEN=default_string_length), PARAMETER    :: description = "[DIPOLE]"
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_result_type), POINTER                      :: results

      NULLIFY (ref_point, results)
      CALL qs_moment_locop(qs_env, magnetic=.FALSE., nmoments=1, reference=use_mom_ref_zero, &
                           ref_point=ref_point, unit_number=-1, vel_reprs=.FALSE., com_nl=.FALSE.)
      CALL get_qs_env(qs_env, results=results)
      IF (.NOT. test_for_result(results=results, description=description)) &
         CALL cp_abort(__LOCATION__, "Failed to compute MaxwellLink dipole source.")
      CALL get_results(results=results, description=description, values=dipole)

   END SUBROUTINE maxwelllink_compute_dipole

! **************************************************************************************************
!> \brief Computes the integrated periodic velocity-gauge current in atomic units.
!> \param qs_env Quickstep environment.
!> \param current Total electronic plus ionic current.
! **************************************************************************************************
   SUBROUTINE maxwelllink_compute_periodic_current(qs_env, current)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: current

      INTEGER                                            :: i, natom, nspin, spin
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      REAL(KIND=dp)                                      :: trace
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, p_im, p_xyz
      TYPE(dbcsr_type), POINTER                          :: tmp_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp

      NULLIFY (dbcsr_dist, dft_control, matrix_s, mos_new, p_im, p_xyz, particle_set, qs_kind_set, &
               rho, row_blk_sizes, rtp, sab_all, sab_orb, tmp_ao)
      CALL get_qs_env(qs_env, &
                      rtp=rtp, &
                      matrix_s=matrix_s, &
                      rho=rho, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      sab_all=sab_all, sab_orb=sab_orb, &
                      dbcsr_dist=dbcsr_dist, nelectron_spin=nelectron_spin)
      IF (rtp%linear_scaling) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET periodic velocity-gauge current requires MO-based RTP.")

      CALL get_rtp(rtp=rtp, mos_new=mos_new)
      nspin = SIZE(mos_new)/2
      natom = SIZE(particle_set, 1)

      ALLOCATE (first_sgf(natom), last_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)
      ALLOCATE (row_blk_sizes(natom))
      CALL dbcsr_convert_offsets_to_sizes(first_sgf, row_blk_sizes, last_sgf)
      DEALLOCATE (first_sgf, last_sgf)

      CALL dbcsr_allocate_matrix_set(p_xyz, 3)
      ALLOCATE (p_xyz(1)%matrix)
      CALL dbcsr_create(matrix=p_xyz(1)%matrix, name="mxl_p_xyz", dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_antisymmetric, row_blk_size=row_blk_sizes, &
                        col_blk_size=row_blk_sizes, mutable_work=.TRUE.)
      CALL cp_dbcsr_alloc_block_from_nbl(p_xyz(1)%matrix, sab_orb)
      CALL dbcsr_set(p_xyz(1)%matrix, 0.0_dp)
      DO i = 2, 3
         ALLOCATE (p_xyz(i)%matrix)
         CALL dbcsr_copy(p_xyz(i)%matrix, p_xyz(1)%matrix, "mxl_p_xyz")
         CALL dbcsr_set(p_xyz(i)%matrix, 0.0_dp)
      END DO
      CALL build_lin_mom_matrix(qs_env, p_xyz)
      DEALLOCATE (row_blk_sizes)

      CALL qs_rho_get(rho, rho_ao_im=p_im)
      NULLIFY (tmp_ao)
      CALL dbcsr_init_p(tmp_ao)
      CALL dbcsr_create(tmp_ao, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry, &
                        name="mxl_current_tmp")
      CALL cp_dbcsr_alloc_block_from_nbl(tmp_ao, sab_all)

      current(:) = 0.0_dp
      DO i = 1, 3
         DO spin = 1, nspin
            CALL dbcsr_set(tmp_ao, 0.0_dp)
            CALL dbcsr_multiply("T", "N", 1.0_dp, p_im(spin)%matrix, p_xyz(i)%matrix, 0.0_dp, tmp_ao)
            CALL dbcsr_trace(tmp_ao, trace)
            current(i) = current(i) - trace + dft_control%rtp_control%vec_pot(i)*REAL(nelectron_spin(spin), dp)
         END DO
      END DO

      IF (.NOT. dft_control%rtp_control%fixed_ions) CALL maxwelllink_add_ionic_current(qs_env, current)

      CALL dbcsr_deallocate_matrix(tmp_ao)
      CALL dbcsr_deallocate_matrix_set(p_xyz)

   END SUBROUTINE maxwelllink_compute_periodic_current

! **************************************************************************************************
!> \brief Adds the ionic part of d mu / dt for Ehrenfest dynamics.
!> \param qs_env Quickstep environment.
!> \param current Current vector to update.
! **************************************************************************************************
   SUBROUTINE maxwelllink_add_ionic_current(qs_env, current)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)         :: current

      INTEGER                                            :: atom_a, iatom, ikind, natom
      INTEGER, DIMENSION(:), POINTER                     :: list
      REAL(KIND=dp)                                      :: zeff
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (atomic_kind_set, list, particle_set, qs_kind_set)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)
      DO ikind = 1, SIZE(qs_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list, natom=natom)
         CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
         DO iatom = 1, natom
            atom_a = list(iatom)
            current(:) = current(:) + zeff*particle_set(atom_a)%v(:)
         END DO
      END DO

   END SUBROUTINE maxwelllink_add_ionic_current

! **************************************************************************************************
!> \brief Builds compact MaxwellLink JSON metadata.
!> \param time Current density time in atomic units.
!> \param energy Total energy in Hartree.
!> \param kinetic_energy Kinetic energy in Hartree.
!> \param potential_energy Field-free potential/electronic energy in Hartree.
!> \param dipole Dipole reported as mux_au/muy_au/muz_au in atomic units.
!> \param dipole_m Current-density dipole reported as mux_m_au/muy_m_au/muz_m_au.
!> \param extra_json JSON string.
! **************************************************************************************************
   SUBROUTINE maxwelllink_build_json(time, energy, kinetic_energy, potential_energy, dipole, dipole_m, extra_json)
      REAL(KIND=dp), INTENT(IN)                          :: time, energy, kinetic_energy, potential_energy
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: dipole, dipole_m
      CHARACTER(LEN=*), INTENT(OUT)                      :: extra_json

      WRITE (extra_json, '(A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,'// &
                         'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A)') &
         '{"time_au":', time, &
         ',"mux_au":', dipole(1), &
         ',"muy_au":', dipole(2), &
         ',"muz_au":', dipole(3), &
         ',"mux_m_au":', dipole_m(1), &
         ',"muy_m_au":', dipole_m(2), &
         ',"muz_m_au":', dipole_m(3), &
         ',"energy_au":', energy, &
         ',"energy_kin_au":', kinetic_energy, &
         ',"energy_pot_au":', potential_energy, &
         '}'

   END SUBROUTINE maxwelllink_build_json

! **************************************************************************************************
!> \brief Builds periodic velocity-gauge MaxwellLink JSON metadata.
!> \param time Current time in atomic units.
!> \param energy Total energy in Hartree.
!> \param kinetic_energy Kinetic energy in Hartree.
!> \param potential_energy Potential/electronic energy in Hartree.
!> \param dipole Effective current-integrated dipole in atomic units.
!> \param dipole_m Midpoint effective dipole in atomic units.
!> \param current Integrated current, i.e. d mu_eff / dt, in atomic units.
!> \param current_density Integrated current divided by cell volume.
!> \param field_midpoint MaxwellLink midpoint electric field.
!> \param vector_potential MaxwellLink physical endpoint vector potential.
!> \param extra_json JSON string.
! **************************************************************************************************
   SUBROUTINE maxwelllink_build_json_periodic(time, energy, kinetic_energy, potential_energy, dipole, &
                                              dipole_m, current, current_density, field_midpoint, &
                                              vector_potential, extra_json)
      REAL(KIND=dp), INTENT(IN)                          :: time, energy, kinetic_energy, potential_energy
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: current, current_density, dipole, dipole_m, &
                                                            field_midpoint, vector_potential
      CHARACTER(LEN=*), INTENT(OUT)                      :: extra_json

      WRITE (extra_json, '(A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,'// &
                          'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,'// &
                          'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,'// &
                          'A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A)') &
         '{"time_au":', time, &
         ',"mux_au":', dipole(1), &
         ',"muy_au":', dipole(2), &
         ',"muz_au":', dipole(3), &
         ',"mux_m_au":', dipole_m(1), &
         ',"muy_m_au":', dipole_m(2), &
         ',"muz_m_au":', dipole_m(3), &
         ',"energy_au":', energy, &
         ',"energy_kin_au":', kinetic_energy, &
         ',"energy_pot_au":', potential_energy, &
         ',"jx_au":', current(1), &
         ',"jy_au":', current(2), &
         ',"jz_au":', current(3), &
         ',"jx_density_au":', current_density(1), &
         ',"jy_density_au":', current_density(2), &
         ',"jz_density_au":', current_density(3), &
         ',"ex_mid_au":', field_midpoint(1), &
         ',"ey_mid_au":', field_midpoint(2), &
         ',"ez_mid_au":', field_midpoint(3), &
         ',"ax_au":', vector_potential(1), &
         ',"ay_au":', vector_potential(2), &
         ',"az_au":', vector_potential(3), &
         ',"periodic_velocity_gauge":true,"energy_field_free_available":false}'

   END SUBROUTINE maxwelllink_build_json_periodic

! **************************************************************************************************
!> \brief Resolves the configured MaxwellLink socket host/path.
!> \param control MaxwellLink socket input control.
!> \return TCP host or UNIX socket path.
! **************************************************************************************************
   FUNCTION maxwelllink_socket_host(control) RESULT(host)
      TYPE(mxl_socket_control_type), INTENT(IN)          :: control
      CHARACTER(LEN=default_path_length)                 :: host

      CHARACTER(LEN=default_path_length)                 :: name

      host = control%host
      IF (.NOT. control%unix) RETURN

      IF (TRIM(control%file) /= "") THEN
         name = control%file
      ELSE
         name = control%host
      END IF

      IF (TRIM(name) /= "" .AND. name(1:1) == "/") THEN
         host = name
      ELSE
         host = TRIM(control%prefix)//TRIM(name)
      END IF

   END FUNCTION maxwelllink_socket_host

END MODULE maxwelllink
