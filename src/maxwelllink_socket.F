!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief MaxwellLink socket protocol helpers.
! **************************************************************************************************
MODULE maxwelllink_socket
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_INT,&
                                              C_NULL_CHAR
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
#ifndef __NO_SOCKETS
   USE sockets_interface,               ONLY: close_socket,&
                                              open_connect_socket,&
                                              readbuffer,&
                                              writebuffer
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'maxwelllink_socket'
   INTEGER, PARAMETER, PRIVATE          :: mxl_msglen = 12

   TYPE maxwelllink_socket_type
      INTEGER                              :: socket = 0
      INTEGER                              :: verbosity = 0
      INTEGER                              :: molecule_id = -1
      REAL(KIND=dp)                        :: init_dt = -1.0_dp
      LOGICAL                              :: connected = .FALSE.
      LOGICAL                              :: handshake = .FALSE.
      CHARACTER(LEN=default_string_length) :: init_json = ""
   END TYPE maxwelllink_socket_type

   PUBLIC :: maxwelllink_socket_type,&
             maxwelllink_socket_connect,&
             maxwelllink_socket_receive_field,&
             maxwelllink_socket_send_source,&
             maxwelllink_socket_shutdown

#ifndef __NO_SOCKETS
   INTERFACE
      SUBROUTINE maxwelllink_open_connect_socket_unix(psockfd, host) &
         BIND(C, NAME="maxwelllink_open_connect_socket_unix")
         IMPORT :: C_CHAR, C_INT
         INTEGER(KIND=C_INT)                      :: psockfd
         CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host
      END SUBROUTINE maxwelllink_open_connect_socket_unix
   END INTERFACE
#endif

CONTAINS

! **************************************************************************************************
!> \brief Connects to a MaxwellLink socket hub.
!> \param mxl_socket MaxwellLink socket state.
!> \param host TCP host name or UNIX socket path.
!> \param port TCP port.
!> \param unix Use UNIX socket if true.
!> \param verbosity Socket logging verbosity.
! **************************************************************************************************
   SUBROUTINE maxwelllink_socket_connect(mxl_socket, host, port, unix, verbosity)
      TYPE(maxwelllink_socket_type), INTENT(INOUT)       :: mxl_socket
      CHARACTER(LEN=*), INTENT(IN)                       :: host
      INTEGER, INTENT(IN)                                :: port
      LOGICAL, INTENT(IN)                                :: unix
      INTEGER, INTENT(IN)                                :: verbosity

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_socket_connect'

#ifdef __NO_SOCKETS
      INTEGER                                            :: handle

      MARK_USED(mxl_socket)
      MARK_USED(host)
      MARK_USED(port)
      MARK_USED(unix)
      MARK_USED(verbosity)
      CALL timeset(routineN, handle)
      CALL cp_abort(__LOCATION__, "MaxwellLink sockets require CP2K to be built with socket support.")
      CALL timestop(handle)
#else
      CHARACTER(LEN=default_path_length)                 :: c_host
      INTEGER                                            :: handle, inet

      CALL timeset(routineN, handle)

      IF (mxl_socket%connected) THEN
         CALL timestop(handle)
         RETURN
      END IF

      mxl_socket%verbosity = verbosity
      c_host = TRIM(host)//C_NULL_CHAR
      IF (unix) THEN
         CALL maxwelllink_open_connect_socket_unix(mxl_socket%socket, c_host)
      ELSE
         inet = 1
         CALL open_connect_socket(mxl_socket%socket, inet, port, c_host)
      END IF
      mxl_socket%connected = .TRUE.

      CALL timestop(handle)
#endif

   END SUBROUTINE maxwelllink_socket_connect

! **************************************************************************************************
!> \brief Receives one MaxwellLink field vector.
!> \param mxl_socket MaxwellLink socket state.
!> \param have_result Whether source data are already available.
!> \param field Uniform electric field in atomic units.
!> \param stop MaxwellLink requested termination.
!> \param received_init INIT packet was received during this call.
!> \param init_dt INIT dt_au value, or negative if absent.
!> \param molecule_id INIT molecule id, or negative if no INIT was received.
! **************************************************************************************************
   SUBROUTINE maxwelllink_socket_receive_field(mxl_socket, have_result, field, stop, &
                                               received_init, init_dt, molecule_id)
      TYPE(maxwelllink_socket_type), INTENT(INOUT)       :: mxl_socket
      LOGICAL, INTENT(IN)                                :: have_result
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: field
      LOGICAL, INTENT(OUT)                               :: stop, received_init
      REAL(KIND=dp), INTENT(OUT)                         :: init_dt
      INTEGER, INTENT(OUT)                               :: molecule_id

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_socket_receive_field'

#ifdef __NO_SOCKETS
      INTEGER                                            :: handle

      MARK_USED(mxl_socket)
      MARK_USED(have_result)
      CALL timeset(routineN, handle)
      field = 0.0_dp
      stop = .TRUE.
      received_init = .FALSE.
      init_dt = -1.0_dp
      molecule_id = -1
      CALL cp_abort(__LOCATION__, "MaxwellLink sockets require CP2K to be built with socket support.")
      CALL timestop(handle)
#else
      CHARACTER(len=mxl_msglen)                          :: header, reply
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (.NOT. mxl_socket%connected) &
         CALL cp_abort(__LOCATION__, "MaxwellLink socket receive requested before connect.")

      field = 0.0_dp
      stop = .FALSE.
      received_init = .FALSE.
      init_dt = -1.0_dp
      molecule_id = -1

      DO
         CALL readbuffer(mxl_socket%socket, header, mxl_msglen)
         SELECT CASE (TRIM(header))
         CASE ("STATUS")
            IF (.NOT. mxl_socket%handshake) THEN
               reply = "NEEDINIT"
            ELSE IF (have_result) THEN
               reply = "HAVEDATA"
            ELSE
               reply = "READY"
            END IF
            CALL writebuffer(mxl_socket%socket, reply, mxl_msglen)
         CASE ("INIT")
            CALL receive_init(mxl_socket)
            received_init = .TRUE.
            init_dt = mxl_socket%init_dt
            molecule_id = mxl_socket%molecule_id
         CASE ("POSDATA", "FIELDDATA")
            CALL receive_field_data(mxl_socket, field)
            EXIT
         CASE ("STOP", "EXIT")
            reply = "BYE"
            CALL writebuffer(mxl_socket%socket, reply, mxl_msglen)
            stop = .TRUE.
            EXIT
         CASE DEFAULT
            CALL cp_abort(__LOCATION__, "Unexpected MaxwellLink message while receiving field: "//TRIM(header))
         END SELECT
      END DO

      CALL timestop(handle)
#endif

   END SUBROUTINE maxwelllink_socket_receive_field

! **************************************************************************************************
!> \brief Sends one MaxwellLink source-current vector.
!> \param mxl_socket MaxwellLink socket state.
!> \param energy Total energy in Hartree.
!> \param source Source-current vector in atomic units.
!> \param extra_json Extra JSON metadata.
!> \param stop MaxwellLink requested termination.
! **************************************************************************************************
   SUBROUTINE maxwelllink_socket_send_source(mxl_socket, energy, source, extra_json, stop)
      TYPE(maxwelllink_socket_type), INTENT(INOUT)       :: mxl_socket
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: source
      CHARACTER(LEN=*), INTENT(IN)                       :: extra_json
      LOGICAL, INTENT(OUT)                               :: stop

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwelllink_socket_send_source'

#ifdef __NO_SOCKETS
      INTEGER                                            :: handle

      MARK_USED(mxl_socket)
      MARK_USED(energy)
      MARK_USED(source)
      MARK_USED(extra_json)
      CALL timeset(routineN, handle)
      stop = .TRUE.
      CALL cp_abort(__LOCATION__, "MaxwellLink sockets require CP2K to be built with socket support.")
      CALL timestop(handle)
#else
      CHARACTER(len=mxl_msglen)                          :: header, reply
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (.NOT. mxl_socket%connected) &
         CALL cp_abort(__LOCATION__, "MaxwellLink socket send requested before connect.")

      stop = .FALSE.
      DO
         CALL readbuffer(mxl_socket%socket, header, mxl_msglen)
         SELECT CASE (TRIM(header))
         CASE ("STATUS")
            reply = "HAVEDATA"
            CALL writebuffer(mxl_socket%socket, reply, mxl_msglen)
         CASE ("GETFORCE", "GETSOURCE")
            CALL send_source_ready(mxl_socket, energy, source, extra_json)
            EXIT
         CASE ("STOP", "EXIT")
            reply = "BYE"
            CALL writebuffer(mxl_socket%socket, reply, mxl_msglen)
            stop = .TRUE.
            EXIT
         CASE DEFAULT
            CALL cp_abort(__LOCATION__, "Unexpected MaxwellLink message while sending source: "//TRIM(header))
         END SELECT
      END DO

      CALL timestop(handle)
#endif

   END SUBROUTINE maxwelllink_socket_send_source

! **************************************************************************************************
!> \brief Closes the MaxwellLink socket.
!> \param mxl_socket MaxwellLink socket state.
! **************************************************************************************************
   SUBROUTINE maxwelllink_socket_shutdown(mxl_socket)
      TYPE(maxwelllink_socket_type), INTENT(INOUT)       :: mxl_socket

#ifndef __NO_SOCKETS
      IF (mxl_socket%connected) CALL close_socket(mxl_socket%socket)
#endif
      mxl_socket%socket = 0
      mxl_socket%connected = .FALSE.
      mxl_socket%handshake = .FALSE.

   END SUBROUTINE maxwelllink_socket_shutdown

#ifndef __NO_SOCKETS
! **************************************************************************************************
!> \brief Receives INIT metadata.
!> \param mxl_socket MaxwellLink socket state.
! **************************************************************************************************
   SUBROUTINE receive_init(mxl_socket)
      TYPE(maxwelllink_socket_type), INTENT(INOUT)       :: mxl_socket

      CHARACTER(LEN=:), ALLOCATABLE                      :: init_json
      INTEGER                                            :: copy_len, nchar
      LOGICAL                                            :: found_dt

      CALL readbuffer(mxl_socket%socket, mxl_socket%molecule_id)
      CALL readbuffer(mxl_socket%socket, nchar)
      IF (nchar < 0) CALL cp_abort(__LOCATION__, "MaxwellLink INIT payload length is negative.")

      IF (nchar > 0) THEN
         ALLOCATE (CHARACTER(LEN=nchar) :: init_json)
         CALL readbuffer(mxl_socket%socket, init_json, nchar)
      ELSE
         init_json = "{}"
      END IF

      mxl_socket%init_json = ""
      copy_len = MIN(LEN(mxl_socket%init_json), LEN(init_json))
      IF (copy_len > 0) mxl_socket%init_json(1:copy_len) = init_json(1:copy_len)
      CALL extract_json_real(init_json, "dt_au", mxl_socket%init_dt, found_dt)
      IF (.NOT. found_dt) mxl_socket%init_dt = -1.0_dp
      mxl_socket%handshake = .TRUE.

   END SUBROUTINE receive_init

! **************************************************************************************************
!> \brief Receives field data after POSDATA/FIELDDATA was consumed.
!> \param mxl_socket MaxwellLink socket state.
!> \param field Uniform electric field in atomic units.
! **************************************************************************************************
   SUBROUTINE receive_field_data(mxl_socket, field)
      TYPE(maxwelllink_socket_type), INTENT(INOUT)       :: mxl_socket
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: field

      INTEGER                                            :: natom
      REAL(KIND=dp), DIMENSION(9)                        :: cell, inv_cell

      CALL readbuffer(mxl_socket%socket, cell, 9)
      CALL readbuffer(mxl_socket%socket, inv_cell, 9)
      CALL readbuffer(mxl_socket%socket, natom)
      IF (natom /= 1) &
         CALL cp_abort(__LOCATION__, "MaxwellLink field packet must contain exactly one vector.")
      CALL readbuffer(mxl_socket%socket, field, 3)

   END SUBROUTINE receive_field_data

! **************************************************************************************************
!> \brief Sends the FORCEREADY/SOURCEREADY payload.
!> \param mxl_socket MaxwellLink socket state.
!> \param energy Total energy in Hartree.
!> \param source Source-current vector in atomic units.
!> \param extra_json Extra JSON metadata.
! **************************************************************************************************
   SUBROUTINE send_source_ready(mxl_socket, energy, source, extra_json)
      TYPE(maxwelllink_socket_type), INTENT(INOUT)       :: mxl_socket
      REAL(KIND=dp), INTENT(IN)                          :: energy
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: source
      CHARACTER(LEN=*), INTENT(IN)                       :: extra_json

      CHARACTER(len=mxl_msglen)                          :: header
      INTEGER                                            :: extra_len, natom
      REAL(KIND=dp), DIMENSION(9)                        :: virial

      header = "FORCEREADY"
      natom = 1
      virial = 0.0_dp
      extra_len = LEN_TRIM(extra_json)

      CALL writebuffer(mxl_socket%socket, header, mxl_msglen)
      CALL writebuffer(mxl_socket%socket, energy)
      CALL writebuffer(mxl_socket%socket, natom)
      CALL writebuffer(mxl_socket%socket, source, 3)
      CALL writebuffer(mxl_socket%socket, virial, 9)
      CALL writebuffer(mxl_socket%socket, extra_len)
      IF (extra_len > 0) CALL writebuffer(mxl_socket%socket, extra_json(1:extra_len), extra_len)

   END SUBROUTINE send_source_ready
#endif

! **************************************************************************************************
!> \brief Extracts a real scalar from a small JSON object without a JSON dependency.
!> \param json JSON text.
!> \param key Key to find.
!> \param value Parsed value.
!> \param found Whether the value was found and parsed.
! **************************************************************************************************
   SUBROUTINE extract_json_real(json, key, value, found)
      CHARACTER(LEN=*), INTENT(IN)                       :: json
      CHARACTER(LEN=*), INTENT(IN)                       :: key
      REAL(KIND=dp), INTENT(OUT)                         :: value
      LOGICAL, INTENT(OUT)                               :: found

      CHARACTER(LEN=:), ALLOCATABLE                      :: token
      INTEGER                                            :: first, last, pos, rel_colon, rel_end, stat

      found = .FALSE.
      value = 0.0_dp

      pos = INDEX(json, '"'//TRIM(key)//'"')
      IF (pos == 0) pos = INDEX(json, TRIM(key))
      IF (pos == 0) RETURN

      rel_colon = INDEX(json(pos:), ":")
      IF (rel_colon == 0) RETURN

      first = pos + rel_colon
      DO WHILE (first <= LEN(json))
         IF (json(first:first) /= " " .AND. json(first:first) /= '"') EXIT
         first = first + 1
      END DO
      IF (first > LEN(json)) RETURN

      rel_end = SCAN(json(first:), ",}")
      IF (rel_end == 0) THEN
         last = LEN(json)
      ELSE
         last = first + rel_end - 2
      END IF
      IF (last < first) RETURN

      token = json(first:last)
      READ (token, *, IOSTAT=stat) value
      found = (stat == 0)

   END SUBROUTINE extract_json_real

END MODULE maxwelllink_socket
