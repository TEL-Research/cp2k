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
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: rtp_control_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_result_methods,               ONLY: get_results,&
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE input_constants,                 ONLY: do_em,&
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
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_moments,                      ONLY: qs_moment_locop
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'maxwelllink'

   TYPE maxwelllink_type
      TYPE(maxwelllink_socket_type)          :: socket
      LOGICAL                                :: active = .FALSE.
      LOGICAL                                :: have_result = .FALSE.
      REAL(KIND=dp), DIMENSION(3)            :: dipole_initial = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: dipole_last = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: field_endpoint = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)            :: field_midpoint = 0.0_dp
      REAL(KIND=dp)                          :: energy_last = 0.0_dp
   END TYPE maxwelllink_type

   PUBLIC :: maxwelllink_type,&
             maxwelllink_init,&
             maxwelllink_receive_field,&
             maxwelllink_send_source,&
             maxwelllink_shutdown,&
             maxwelllink_validate_rtp_control

CONTAINS

! **************************************************************************************************
!> \brief Validates RT-TDDFT control settings for MaxwellLink coupling.
!> \param rtp_control RTP control settings.
! **************************************************************************************************
   SUBROUTINE maxwelllink_validate_rtp_control(rtp_control)
      TYPE(rtp_control_type), INTENT(IN)                 :: rtp_control

      IF (.NOT. rtp_control%mxl_socket%active) RETURN

      IF (.NOT. rtp_control%fixed_ions) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET is currently implemented only for fixed-ion RT-TDDFT.")
      IF (rtp_control%velocity_gauge) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET currently supports only non-periodic length-gauge RT-TDDFT.")
      IF (rtp_control%propagator /= do_em) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET requires REAL_TIME_PROPAGATION%PROPAGATOR EM.")

   END SUBROUTINE maxwelllink_validate_rtp_control

! **************************************************************************************************
!> \brief Initializes MaxwellLink coupling.
!> \param mxl MaxwellLink runtime state.
!> \param qs_env Quickstep environment.
!> \param rtp_control RTP control settings.
!> \param energy Initial total energy.
! **************************************************************************************************
   SUBROUTINE maxwelllink_init(mxl, qs_env, rtp_control, energy)
      TYPE(maxwelllink_type), INTENT(INOUT)              :: mxl
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rtp_control_type), INTENT(INOUT)              :: rtp_control
      REAL(KIND=dp), INTENT(IN)                          :: energy

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_init'

      CHARACTER(LEN=default_path_length)                 :: socket_host
      INTEGER                                            :: handle, output_unit
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      IF (.NOT. rtp_control%mxl_socket%active) THEN
         CALL timestop(handle)
         RETURN
      END IF

      NULLIFY (cell, para_env)
      CALL get_qs_env(qs_env, cell=cell, para_env=para_env)
      IF (.NOT. rtp_control%fixed_ions) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET is currently implemented only for fixed-ion RT-TDDFT.")
      IF (rtp_control%velocity_gauge) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET currently supports only non-periodic length-gauge RT-TDDFT.")
      IF (ANY(cell%perd(1:3) /= 0)) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET currently supports only non-periodic cells.")
      IF (rtp_control%propagator /= do_em) &
         CALL cp_abort(__LOCATION__, "MXL_SOCKET requires REAL_TIME_PROPAGATION%PROPAGATOR EM.")

      mxl%active = .TRUE.
      mxl%have_result = .FALSE.
      mxl%energy_last = energy
      mxl%field_endpoint = 0.0_dp
      mxl%field_midpoint = 0.0_dp
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

      IF (rtp_control%mxl_socket%reset_dipole) THEN
         CALL maxwelllink_compute_dipole(qs_env, mxl%dipole_initial)
      ELSE
         mxl%dipole_initial = 0.0_dp
      END IF
      CALL maxwelllink_compute_dipole(qs_env, mxl%dipole_last)
      mxl%dipole_last = mxl%dipole_last - mxl%dipole_initial

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
         ! MaxwellLink sends midpoint fields; EM will average this endpoint with the previous one.
         mxl%field_midpoint = field
         field = 2.0_dp*mxl%field_midpoint - mxl%field_endpoint
         mxl%field_endpoint = field
         rtp_control%field = field
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
!> \param energy Current total energy.
!> \param stop MaxwellLink requested termination.
! **************************************************************************************************
   SUBROUTINE maxwelllink_send_source(mxl, qs_env, dt, sim_time, energy, stop)
      TYPE(maxwelllink_type), INTENT(INOUT)              :: mxl
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: dt, sim_time, energy
      LOGICAL, INTENT(OUT)                               :: stop

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_send_source'

      CHARACTER(LEN=max_line_length)                     :: extra_json
      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: energy_field_free, source_time
      REAL(KIND=dp), DIMENSION(3)                        :: dipole_abs, dipole_end, dipole_json, source
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      stop = .FALSE.
      IF (.NOT. mxl%active) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL maxwelllink_compute_dipole(qs_env, dipole_abs)
      dipole_end = dipole_abs - mxl%dipole_initial
      source(:) = (dipole_end(:) - mxl%dipole_last(:))/dt
      dipole_json(:) = dipole_end(:)
      source_time = sim_time
      energy_field_free = energy + DOT_PRODUCT(dipole_abs, mxl%field_endpoint)
      CALL maxwelllink_build_json(source_time, energy_field_free, dipole_json, dipole_end, extra_json)

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)
      IF (para_env%is_source()) THEN
         CALL maxwelllink_socket_send_source(mxl%socket, energy_field_free, source, TRIM(extra_json), stop)
      ELSE
         stop = .FALSE.
      END IF
      CALL para_env%bcast(stop)

      mxl%have_result = .FALSE.
      mxl%dipole_last = dipole_end
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
!> \brief Builds compact MaxwellLink JSON metadata.
!> \param time Current density time in atomic units.
!> \param energy Total energy in Hartree.
!> \param dipole Dipole reported as mux_au/muy_au/muz_au in atomic units.
!> \param dipole_m Current-density dipole reported as mux_m_au/muy_m_au/muz_m_au.
!> \param extra_json JSON string.
! **************************************************************************************************
   SUBROUTINE maxwelllink_build_json(time, energy, dipole, dipole_m, extra_json)
      REAL(KIND=dp), INTENT(IN)                          :: time, energy
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: dipole, dipole_m
      CHARACTER(LEN=*), INTENT(OUT)                      :: extra_json

      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

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
         ',"energy_kin_au":', zero, &
         ',"energy_pot_au":', energy, &
         '}'

   END SUBROUTINE maxwelllink_build_json

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
