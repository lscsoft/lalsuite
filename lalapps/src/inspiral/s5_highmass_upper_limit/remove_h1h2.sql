DELETE FROM
	coinc_inspiral
WHERE
	ifos == "H1,H2";
	
DELETE FROM
	coinc_event
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_inspiral
	);

DELETE FROM
	coinc_event_map
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
	);

DELETE FROM
	sngl_inspiral
WHERE
	event_id NOT IN (
		SELECT
			event_id
		FROM
			coinc_event_map
		WHERE
			table_name == "sngl_inspiral"
	);
