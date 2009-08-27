-- remove coincs when H1+H2 are the only instruments on

DELETE FROM
	coinc_event
WHERE
	instruments == "H1,H2";

DELETE FROM
	coinc_inspiral
WHERE
	ifos == "H1,H2";

-- remove H2+L1 coincs when H1+H2+L1 are on

DELETE FROM
	coinc_inspiral
WHERE
	coinc_event_id IN (
		SELECT
			coinc_inspiral.coinc_event_id
		FROM
			coinc_inspiral
			JOIN coinc_event ON (
				coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
			)
		WHERE
			coinc_inspiral.ifos == "H2,L1"
			AND coinc_event.instruments == "H1,H2,L1"
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
