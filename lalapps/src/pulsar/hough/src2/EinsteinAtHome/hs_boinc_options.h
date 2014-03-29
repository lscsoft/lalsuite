/*
*  Copyright (C) 2007 Bernd Machenschalk
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* C++ -> C interface for BOINC_OPTION handling etc.
*/

#define EAH_SHMEM_APP_NAME "EinsteinHS"
#define EAH_SHMEM_SIZE 1024

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double skypos_rac;
  double skypos_dec;
  double frequency;
  double bandwidth;
  double cand_frequency;
  double cand_spindown;
  double cand_rac;
  double cand_dec;
  double cand_hough_sign;
} t_progress;

extern t_progress boincv6_progress;

#ifdef BOINC_APIV6
extern int setup_shmem(void);
#endif

extern void set_boinc_options(void);

extern int eah_userid, eah_hostid;
extern char*eah_hostcpid, *eah_username;

#ifdef __cplusplus
}
#endif
