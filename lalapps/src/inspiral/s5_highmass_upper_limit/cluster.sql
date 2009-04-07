SELECT
	"Number of coincs before clustering: " || count(*)
FROM
	coinc_event;

CREATE TEMPORARY TABLE maxsnrs AS
	SELECT
		coinc_inspiral.end_time AS end_time,
		coinc_event.time_slide_id AS time_slide_id,
		MAX(coinc_inspiral.snr) AS snr
	FROM
		coinc_inspiral
		JOIN coinc_event ON (
			coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
		)
	GROUP BY
		coinc_inspiral.end_time,
		coinc_event.time_slide_id;

CREATE INDEX ms_et_index ON maxsnrs (end_time, time_slide_id);
DELETE FROM
	coinc_inspiral
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_inspiral.coinc_event_id
		FROM
			coinc_inspiral
			JOIN coinc_event ON (
				coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
			)
			JOIN maxsnrs ON (
				maxsnrs.end_time == coinc_inspiral.end_time
				AND maxsnrs.time_slide_id == coinc_event.time_slide_id
			)
		WHERE
			coinc_inspiral.snr >= maxsnrs.snr
	);
DROP INDEX ms_et_index;
DROP TABLE maxsnrs;

DELETE FROM
	coinc_event
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_inspiral
	);

-- find coincs that are within 10 s of coincs of higher effective SNR that involved the same instruments and were found in the same time slide
CREATE INDEX ci_et_index ON coinc_inspiral (end_time);
DELETE FROM
	coinc_event
WHERE
	coinc_event_id IN (
		SELECT DISTINCT
			coinc_inspiral_a.coinc_event_id
		FROM
			coinc_inspiral AS coinc_inspiral_a
			JOIN coinc_inspiral AS coinc_inspiral_b INDEXED BY ci_et_index ON (
				coinc_inspiral_b.end_time BETWEEN coinc_inspiral_a.end_time - 10 AND coinc_inspiral_a.end_time + 10
				AND coinc_inspiral_b.coinc_event_id != coinc_inspiral_a.coinc_event_id
			)
			JOIN coinc_event AS coinc_event_a INDEXED BY sqlite_autoindex_coinc_event_1 ON (
				coinc_event_a.coinc_event_id == coinc_inspiral_a.coinc_event_id
			)
			JOIN coinc_event AS coinc_event_b INDEXED BY sqlite_autoindex_coinc_event_1 ON (
				coinc_event_b.coinc_event_id == coinc_inspiral_b.coinc_event_id
			)
		WHERE
			coinc_event_b.time_slide_id == coinc_event_a.time_slide_id
			AND coinc_inspiral_b.snr >= coinc_inspiral_a.snr
			AND coinc_event_b.instruments == coinc_event_a.instruments
			AND coinc_event_b.coinc_def_id == coinc_event_a.coinc_def_id
			AND abs((coinc_inspiral_b.end_time - coinc_inspiral_a.end_time) + (coinc_inspiral_b.end_time_ns - coinc_inspiral_a.end_time_ns) * 1e-9) < 10.0
	);
DROP INDEX ci_et_index;

SELECT
	"Number of coincs after clustering: " || count(*)
FROM
	coinc_event;

-- delete unused coinc_inspiral rows
DELETE FROM
	coinc_inspiral
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
	);

-- delete unused coinc_event_map rows
DELETE FROM
	coinc_event_map
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
	);

-- delete unused sngl_inspiral rows
DELETE FROM
	sngl_inspiral
WHERE
	event_id NOT IN (
		SELECT
			event_id
		FROM
			coinc_event_map
		WHERE
			table_name == 'sngl_inspiral'
	);
