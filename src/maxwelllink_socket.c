/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef __NO_SOCKETS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/un.h>

void maxwelllink_open_connect_socket_unix(int *psockfd, char *host) {
  int sockfd;
  struct sockaddr_un serv_addr;

  memset(&serv_addr, 0, sizeof(serv_addr));
  serv_addr.sun_family = AF_UNIX;
  if (strlen(host) >= sizeof(serv_addr.sun_path)) {
    fprintf(stderr, "Error opening MaxwellLink UNIX socket: path too long\n");
    exit(-1);
  }
  strcpy(serv_addr.sun_path, host);

  sockfd = socket(AF_UNIX, SOCK_STREAM, 0);
  if (sockfd < 0) {
    perror("Error opening MaxwellLink UNIX socket");
    exit(-1);
  }

  if (connect(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
    perror("Error opening MaxwellLink UNIX socket");
    exit(-1);
  }

  *psockfd = sockfd;
}

#endif
