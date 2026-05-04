!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief MaxwellLink control data types.
! **************************************************************************************************
MODULE maxwelllink_types
   USE kinds,                           ONLY: default_path_length
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'maxwelllink_types'

   INTEGER, PARAMETER                    :: mxl_periodic_dipole_berry = 1
   INTEGER, PARAMETER                    :: mxl_periodic_dipole_current = 2

   TYPE mxl_socket_control_type
      LOGICAL                              :: active = .FALSE.
      LOGICAL                              :: unix = .FALSE.
      LOGICAL                              :: reset_dipole = .FALSE.
      INTEGER                              :: molecule_id = -1
      INTEGER                              :: periodic_dipole_method = mxl_periodic_dipole_berry
      INTEGER                              :: port = 31415
      INTEGER                              :: verbosity = 0
      CHARACTER(LEN=default_path_length)   :: host = "localhost"
      CHARACTER(LEN=default_path_length)   :: file = ""
      CHARACTER(LEN=default_path_length)   :: prefix = "/tmp/socketmxl_"
   END TYPE mxl_socket_control_type

   PUBLIC :: mxl_periodic_dipole_berry,&
             mxl_periodic_dipole_current,&
             mxl_socket_control_type

END MODULE maxwelllink_types
