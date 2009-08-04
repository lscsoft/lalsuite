--
-- coinc_definer clean up
--


CREATE TABLE _idmap_ AS
	SELECT
		old_definer.coinc_def_id AS old,
		MIN(new_definer.coinc_def_id) AS new
	FROM
		coinc_definer AS old_definer
		JOIN coinc_definer AS new_definer ON (
			new_definer.search == old_definer.search
			AND new_definer.search_coinc_type == old_definer.search_coinc_type
		)
	GROUP BY
		old_definer.coinc_def_id;

CREATE INDEX idm_o_index ON _idmap_ (old);
UPDATE coinc_event SET coinc_def_id = (SELECT new FROM _idmap_ WHERE old == coinc_def_id);
DROP INDEX idm_o_index;

DELETE FROM coinc_definer WHERE coinc_def_id IN (SELECT old FROM _idmap_ WHERE old != new);

DROP TABLE _idmap_;


--
-- time_slide clean up
--


CREATE TEMPORARY TABLE time_slide_ids AS SELECT DISTINCT time_slide_id AS time_slide_id FROM time_slide;
CREATE INDEX ts_ioid_index ON time_slide (instrument, offset, time_slide_id);
CREATE INDEX ts_io_index ON time_slide (instrument, offset);
CREATE INDEX ts_iid_index ON time_slide (instrument, time_slide_id);
CREATE INDEX ts_oid_index ON time_slide (offset, time_slide_id);
CREATE INDEX ts_i_index ON time_slide (instrument);
CREATE INDEX ts_o_index ON time_slide (offset);
CREATE INDEX ts_o_index ON time_slide (time_slide_id);
CREATE INDEX tsid_id_index ON time_slide_ids (time_slide_id);

CREATE TABLE _idmap_ AS
	SELECT
		old_ids.time_slide_id AS old,
		(
			SELECT
				MIN(new_ids.time_slide_id)
			FROM
				time_slide_ids AS new_ids
			WHERE
				NOT EXISTS (
					SELECT
						*
					FROM
						time_slide AS alink
						JOIN time_slide AS blink
					WHERE
						alink.time_slide_id == old_ids.time_slide_id
						AND blink.time_slide_id == new_ids.time_slide_id
						AND blink.instrument == alink.instrument
						AND blink.offset != alink.offset
				)
		) AS new
	FROM
		time_slide_ids AS old_ids;

DROP INDEX ts_ioid_index;
DROP INDEX ts_io_index;
DROP INDEX ts_iid_index;
DROP INDEX ts_oid_index;
DROP INDEX ts_i_index;
DROP INDEX ts_o_index;
DROP INDEX ts_id_index;
DROP INDEX tsid_id_index;


DROP TABLE time_slide_ids;

CREATE INDEX idm_o_index ON _idmap_ (old);
UPDATE coinc_event SET time_slide_id = (SELECT new FROM _idmap_ WHERE old == time_slide_id);
DROP INDEX idm_o_index;

DELETE FROM time_slide WHERE time_slide_id IN (SELECT old FROM _idmap_ WHERE old != new);

DROP TABLE _idmap_;
