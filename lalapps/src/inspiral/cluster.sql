PRAGMA journal_mode = MEMORY;
PRAGMA locking_mode = EXCLUSIVE;
PRAGMA synchronous = OFF;
PRAGMA temp_store = MEMORY;
-- PRAGMA temp_store_directory = '/tmp';

SELECT
	"Number of coincs before clustering: " || count(*)
FROM
	coinc_event;

--
-- construct a look-up table of playground/non-playground state
--

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
CREATE INDEX tmpindex1 ON is_playground (coinc_event_id);

--
-- create a look-up table of info required for clustering
--

CREATE TEMPORARY TABLE _cluster_info_ AS
	SELECT
		coinc_event.coinc_event_id AS coinc_event_id,
		(coinc_event.time_slide_id || ";" || coinc_event.instruments || ";" || is_playground.is_playground) AS category,
		(coinc_inspiral.end_time - (SELECT MIN(end_time) FROM coinc_inspiral)) + 1e-9 * coinc_inspiral.end_time_ns AS end_time,
		coinc_inspiral.snr AS snr
	FROM
		coinc_event
		JOIN is_playground ON (
			is_playground.coinc_event_id == coinc_event.coinc_event_id
		)
		JOIN coinc_inspiral ON (
			coinc_inspiral.coinc_event_id == coinc_event.coinc_event_id
		);
DROP INDEX tmpindex1;
DROP TABLE is_playground;
CREATE INDEX tmpindex1 ON _cluster_info_ (coinc_event_id);
CREATE INDEX tmpindex2 ON _cluster_info_ (category, end_time, snr);

--
-- delete coincs that are within 10 s of coincs with higher SNR in the same
-- category
--

DELETE FROM
	coinc_event
WHERE
	EXISTS (
		SELECT
			*
		FROM
			_cluster_info_ AS _cluster_info_a_
			JOIN _cluster_info_ AS _cluster_info_b_ ON (
				_cluster_info_b_.category == _cluster_info_a_.category
				AND (_cluster_info_b_.end_time BETWEEN _cluster_info_a_.end_time - 10.0 AND _cluster_info_a_.end_time + 10.0)
				AND _cluster_info_b_.snr > _cluster_info_a_.snr
			)
		WHERE
			_cluster_info_a_.coinc_event_id == coinc_event.coinc_event_id
	);
DROP INDEX tmpindex1;
DROP INDEX tmpindex2;
DROP TABLE _cluster_info_;

SELECT
	"Number of coincs after clustering: " || count(*)
FROM
	coinc_event;

--
-- delete unused coinc_inspiral rows
--

DELETE FROM
	coinc_inspiral
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
	);

--
-- delete unused coinc_event_map rows
--

DELETE FROM
	coinc_event_map
WHERE
	coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event
	);

--
-- delete unused sngl_inspiral rows
--

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

--
-- shrink the file
--

-- VACUUM;
