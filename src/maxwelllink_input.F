!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief MaxwellLink input helpers.
! **************************************************************************************************
MODULE maxwelllink_input
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_val_types,                 ONLY: lchar_t
   USE maxwelllink_types,               ONLY: mxl_socket_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'maxwelllink_input'

   PUBLIC :: create_mxl_socket_section,&
             read_mxl_socket_control

CONTAINS

! **************************************************************************************************
!> \brief Creates the MaxwellLink socket subsection.
!> \param section The created MXL_SOCKET section.
! **************************************************************************************************
   SUBROUTINE create_mxl_socket_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MXL_SOCKET", &
                          description="Couples time-dependent dynamics to MaxwellLink through a socket. "// &
                          "For real-time TDDFT this receives a uniform electric field in atomic units and "// &
                          "returns the molecular source current d(mu)/dt in atomic units.", &
                          n_keywords=8, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="HOST", &
                          description="Host name for TCP sockets. If UNIX is true and FILE is empty, this is "// &
                          "interpreted as the UNIX socket path.", &
                          usage="HOST localhost", default_c_val="localhost")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PORT", &
                          description="TCP port for the MaxwellLink socket hub.", &
                          usage="PORT 31415", default_i_val=31415)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNIX", &
                          description="Use a UNIX domain socket instead of a TCP socket.", &
                          usage="UNIX", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FILE", &
                          description="UNIX socket filename or absolute path. Relative names are appended to PREFIX.", &
                          usage="FILE molecule.sock", type_of_var=lchar_t, default_lc_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PREFIX", &
                          description="Prefix used when FILE is a relative UNIX socket name.", &
                          usage="PREFIX /tmp/socketmxl_", type_of_var=lchar_t, &
                          default_lc_val="/tmp/socketmxl_")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOLECULE_ID", &
                          description="Expected MaxwellLink molecule id. Negative values disable checking.", &
                          usage="MOLECULE_ID 0", default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESET_DIPOLE", &
                          description="Subtract the initial dipole from diagnostic dipoles returned to MaxwellLink. "// &
                          "The returned source current is unchanged apart from roundoff.", &
                          usage="RESET_DIPOLE", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VERBOSITY", &
                          description="Verbosity level for MaxwellLink socket logging.", &
                          usage="VERBOSITY 0", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mxl_socket_section

! **************************************************************************************************
!> \brief Reads the MaxwellLink socket subsection into control data.
!> \param rtp_section REAL_TIME_PROPAGATION input section.
!> \param control MaxwellLink socket control data.
! **************************************************************************************************
   SUBROUTINE read_mxl_socket_control(rtp_section, control)
      TYPE(section_vals_type), POINTER                   :: rtp_section
      TYPE(mxl_socket_control_type), INTENT(INOUT)       :: control

      LOGICAL                                            :: is_present
      TYPE(section_vals_type), POINTER                   :: subsection

      subsection => section_vals_get_subs_vals(rtp_section, "MXL_SOCKET")
      CALL section_vals_get(subsection, explicit=is_present)
      control%active = is_present
      IF (.NOT. is_present) RETURN

      CALL section_vals_val_get(subsection, "HOST", c_val=control%host)
      CALL section_vals_val_get(subsection, "PORT", i_val=control%port)
      CALL section_vals_val_get(subsection, "UNIX", l_val=control%unix)
      CALL section_vals_val_get(subsection, "FILE", c_val=control%file)
      CALL section_vals_val_get(subsection, "PREFIX", c_val=control%prefix)
      CALL section_vals_val_get(subsection, "MOLECULE_ID", i_val=control%molecule_id)
      CALL section_vals_val_get(subsection, "RESET_DIPOLE", l_val=control%reset_dipole)
      CALL section_vals_val_get(subsection, "VERBOSITY", i_val=control%verbosity)

   END SUBROUTINE read_mxl_socket_control

END MODULE maxwelllink_input
