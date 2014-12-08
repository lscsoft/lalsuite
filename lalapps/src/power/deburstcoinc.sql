-- Copyright (C) 2014 Kipp Cannon
--
-- This program is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the
-- Free Software Foundation; either version 2 of the License, or (at your
-- option) any later version.
--
-- This program is distributed in the hope that it will be useful, but
-- WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
-- Public License for more details.
--
-- You should have received a copy of the GNU General Public License along
-- with this program; if not, write to the Free Software Foundation, Inc.,
-- 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

-- SQL script to remove burst-related coinc tables
-- before running this, run debinjfind.sql and deburca.sql to clean up
-- process metadata

DELETE TABLE coinc_event;
DELETE TABLE coinc_event_map;
DELETE TABLE coinc_definer;
DELETE TABLE multi_burst;
