 /*
  * Copyright (C) 2004, 2005 Cristina V. Torres
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
  *
  */

 /*
  * Author: Cristina Valeria Torres (LIGO Livingston)
  */

#ifndef FRAMEDATACONVERTERREAL4_H
#define FRAMEDATACONVERTERREAL4_H

#include "tracksearch.h"
#include "tracksearchAverager.h"
#include "tracksearchToolbox.h"
#include <unistd.h>
#include <errno.h>

define(`TYPECODE',`D')
include(`CreateConverterTypeDef.m4')

define(`TYPECODE',`S')
include(`CreateConverterTypeDef.m4')

define(`TYPECODE',`I2')
include(`CreateConverterTypeDef.m4')

define(`TYPECODE',`I4')
include(`CreateConverterTypeDef.m4')

define(`TYPECODE',`I8')
include(`CreateConverterTypeDef.m4')

define(`TYPECODE',`U2')
include(`CreateConverterTypeDef.m4')

define(`TYPECODE',`U4')
include(`CreateConverterTypeDef.m4')

define(`TYPECODE',`U8')
include(`CreateConverterTypeDef.m4')

#endif
