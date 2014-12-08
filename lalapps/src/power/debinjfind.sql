-- SQL script to remove traces of ligolw_binjfind from a database.
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
			program == 'ligolw_binjfind'
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
			program == 'ligolw_binjfind'
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
					program == 'ligolw_binjfind'
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
			program == 'ligolw_binjfind'
	);
