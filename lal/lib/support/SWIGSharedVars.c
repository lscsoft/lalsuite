//
// Copyright (C) 2022 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

///
/// \defgroup SWIGSharedVars_c File SWIGSharedVars.c
/// \ingroup lal_swig
/// \brief Flag variables shared between the LAL SWIG wrappers
/// \author Karl Wette
///
/// This file contains definitions of internal variables which are shared between the LAL SWIG wrappers.
/// For simple C variables (i.e. not scripting language objects) it is easier to share them as variables
/// in the LALSupport library, which is then dynamically linked to each SWIG wrapper library.
///
/// @{

///
/// The \c swig_lal_do_redirect_stdouterr variable turns on standard output/error redirection for
/// all LAL libraries. The \c swig_lal_has_stdouterr_been_redirected variable indicates where
/// standard standard output/error redirection is currently in force. See <tt>SWIGCommon.i</tt> for
/// further information.
///
int swig_lal_do_redirect_stdouterr = 0;
int swig_lal_has_stdouterr_been_redirected = 0;

/// @}
