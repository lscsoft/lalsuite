PRAGMA temp_store_directory = '/tmp';

SELECT
	"Number of coincs before clustering: " || count(*)
FROM
	coinc_event;

-- construct "is playground" look-up table.  See
-- glue.segmentsUtils.S2playground() for definition of playground segments
CREATE TEMPORARY TABLE is_playground AS
	SELECT
		coinc_inspiral.coinc_event_id AS coinc_event_id,
		-- is this a zero-lag coinc, and did the last playground
		-- segment start less than 600 seconds prior to it?
		NOT EXISTS (SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0) AND ((coinc_inspiral.end_time - 729273613) % 6370) < 600 AS is_playground
	FROM
		coinc_inspiral
		JOIN coinc_event ON (
			coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
		);
CREATE INDEX ip_cei_index ON is_playground (coinc_event_id);


-- temporary table containing highest SNR in each GPS second in each time
-- slide in each on-instruments category for both playground and non-playground
CREATE TEMPORARY TABLE maxsnrs AS
	SELECT
		coinc_inspiral.end_time AS end_time,
		coinc_event.time_slide_id AS time_slide_id,
		coinc_event.instruments AS on_instruments,
		is_playground.is_playground AS is_playground,
		MAX(coinc_inspiral.snr) AS snr
	FROM
		coinc_inspiral
		JOIN coinc_event ON (
			coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
		)
		JOIN is_playground ON (
			is_playground.coinc_event_id == coinc_inspiral.coinc_event_id
		)
	GROUP BY
		coinc_inspiral.end_time,
		coinc_event.time_slide_id,
		coinc_event.instruments,
		is_playground.is_playground;

-- remove coincs whose SNRs are less than the highest SNR in their GPS
-- second and time slide and on-instrument category
CREATE INDEX ms_et_index ON maxsnrs (end_time, time_slide_id, on_instruments, is_playground);
DELETE FROM
	coinc_inspiral
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_inspiral.coinc_event_id
		FROM
			coinc_inspiral
			JOIN coinc_event INDEXED BY sqlite_autoindex_coinc_event_1 ON (
				coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
			)
			JOIN is_playground INDEXED BY ip_cei_index ON (
				is_playground.coinc_event_id == coinc_inspiral.coinc_event_id
			)
			JOIN maxsnrs INDEXED BY ms_et_index ON (
				maxsnrs.end_time == coinc_inspiral.end_time
				AND maxsnrs.time_slide_id == coinc_event.time_slide_id
				AND maxsnrs.on_instruments == coinc_event.instruments
				AND maxsnrs.is_playground == is_playground.is_playground
			)
		WHERE
			coinc_inspiral.snr >= maxsnrs.snr
	);
DROP INDEX ms_et_index;
DROP TABLE maxsnrs;

-- remove unused coinc_event rows
DELETE FROM
	coinc_event
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_inspiral
	);

-- find coincs that are within 10 s of coincs of higher effective SNR that
-- involved the same instruments and were found in the same time slide
-- NOTE:  this code involves some *very* careful use of indeces to make
-- this go quick, including the application of the constraints in a special
-- order
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
				-- apply coarse \Delta t cut
				coinc_inspiral_b.end_time BETWEEN coinc_inspiral_a.end_time - 10 AND coinc_inspiral_a.end_time + 10
				-- don't cluster coincs with themselves
				AND coinc_inspiral_b.coinc_event_id != coinc_inspiral_a.coinc_event_id
			)
			JOIN coinc_event AS coinc_event_a INDEXED BY sqlite_autoindex_coinc_event_1 ON (
				coinc_event_a.coinc_event_id == coinc_inspiral_a.coinc_event_id
			)
			JOIN is_playground AS is_playground_a INDEXED BY ip_cei_index ON (
				is_playground_a.coinc_event_id == coinc_inspiral_a.coinc_event_id
			)
			JOIN coinc_event AS coinc_event_b INDEXED BY sqlite_autoindex_coinc_event_1 ON (
				coinc_event_b.coinc_event_id == coinc_inspiral_b.coinc_event_id
			)
			JOIN is_playground AS is_playground_b INDEXED BY ip_cei_index ON (
				is_playground_b.coinc_event_id == coinc_inspiral_b.coinc_event_id
			)
		WHERE
			coinc_event_b.time_slide_id == coinc_event_a.time_slide_id
			AND coinc_inspiral_b.snr >= coinc_inspiral_a.snr
			AND coinc_event_b.instruments == coinc_event_a.instruments
			AND is_playground_b.is_playground == is_playground_a.is_playground
			AND coinc_event_b.coinc_def_id == coinc_event_a.coinc_def_id
			-- apply fine \Delta t cut
			AND abs((coinc_inspiral_b.end_time - coinc_inspiral_a.end_time) + (coinc_inspiral_b.end_time_ns - coinc_inspiral_a.end_time_ns) * 1e-9) < 10.0
	);
DROP INDEX ci_et_index;
DROP INDEX ip_cei_index;
DROP TABLE is_playground;

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

-- VACUUM, this reduces the file size by ~6
VACUUM;
