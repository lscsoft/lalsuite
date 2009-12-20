PRAGMA journal_mode = MEMORY;
PRAGMA locking_mode = EXCLUSIVE;
PRAGMA synchronous = OFF;
PRAGMA temp_store = MEMORY;
-- PRAGMA temp_store_directory = '/tmp';

-- remove coincs when H1+H2 are the only instruments on

DELETE FROM
	coinc_event
WHERE
	instruments == "H1,H2";

-- remove coincs when H1+H2 are the only participants or when H2+L1 are the
-- only participants when H1+H2+L1 are on

DELETE FROM
	coinc_inspiral
WHERE
	ifos == "H1,H2"
	OR (
		ifos == "H2,L1"
		AND coinc_event_id IN (
			SELECT
				coinc_event_id
			FROM
				coinc_event
			WHERE
				instruments == "H1,H2,L1"
		)
	);

-- remove unused rows from the coinc_event and coinc_inspiral tables

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
	coinc_inspiral
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
	);

-- remove unused rows from the coinc_event_map table

DELETE FROM
	coinc_event_map
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
	);
