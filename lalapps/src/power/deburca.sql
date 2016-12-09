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


-- SQL script to remove traces of lalapps_burca from a database.
-- apply with sqlite3 command interpreter, or lalapps_runsqlite


DELETE FROM
	process_params
WHERE
	process_id IN (
		SELECT
			process_id
		FROM
			process
		WHERE
			program == 'lalapps_burca'
	);
	
DELETE FROM
	search_summary
WHERE
	process_id IN (
		SELECT
			process_id
		FROM
			process
		WHERE
			program == 'lalapps_burca'
	);

DELETE FROM
	coinc_event_map
WHERE
	coinc_event_id IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
		WHERE
			process_id IN (
				SELECT
					process_id
				FROM
					process
				WHERE
					program == 'lalapps_burca'
			)
	);
	
DELETE FROM
	coinc_event
WHERE
	process_id IN (
		SELECT
			process_id
		FROM
			process
		WHERE
			program == 'lalapps_burca'
	);
