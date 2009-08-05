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
CREATE INDEX ce_cdid_index ON coinc_event (coinc_def_id);
UPDATE coinc_event SET coinc_def_id = (SELECT new FROM _idmap_ WHERE old == coinc_def_id);
DROP INDEX idm_o_index;
DROP INDEX ce_cdid_index;

DELETE FROM coinc_definer WHERE coinc_def_id IN (SELECT old FROM _idmap_ WHERE old != new);

DROP TABLE _idmap_;

--
-- time_slide clean up
--

CREATE TEMPORARY TABLE time_slide_ids AS SELECT DISTINCT time_slide_id AS time_slide_id FROM time_slide;
CREATE INDEX tsid_id_index ON time_slide_ids (time_slide_id);

-- BEGIN CHADS CHANGES TO SPEED UP TIME SLIDE ID FIXING
CREATE TABLE new_slide AS
        SELECT time_slide_id, instrument + offset AS vec FROM time_slide;
CREATE INDEX ns_tsid_index ON new_slide (time_slide_id);
CREATE INDEX ns_tsidv_index ON new_slide (time_slide_id, vec);
CREATE INDEX ns_v_index ON new_slide (vec);
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
                                                new_slide AS alink
                                                JOIN new_slide AS blink
                                        WHERE
                                                alink.time_slide_id == old_ids.time_slide_id
                                                AND blink.vec == new_ids.vec
                                )
                ) AS new
        FROM
                time_slide_ids AS old_ids;
DROP INDEX ns_tsid_index;
DROP INDEX ns_tsidv_index;
DROP INDEX ns_v_index;
--- END CHADS CHANGES

DROP TABLE time_slide_ids;

CREATE INDEX idm_o_index ON _idmap_ (old);
UPDATE coinc_event SET time_slide_id = (SELECT new FROM _idmap_ WHERE old == time_slide_id);
DROP INDEX idm_o_index;

DELETE FROM time_slide WHERE time_slide_id IN (SELECT old FROM _idmap_ WHERE old != new);

DROP TABLE _idmap_;
