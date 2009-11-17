PRAGMA journal_mode = MEMORY;
PRAGMA locking_mode = EXCLUSIVE;
PRAGMA synchronous = OFF;
PRAGMA temp_store = MEMORY;
-- PRAGMA temp_store_directory = '/tmp';

--
-- coinc_definer clean up
--

CREATE TEMPORARY TABLE _idmap_ AS
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
CREATE INDEX tmpindex ON _idmap_ (old);

UPDATE coinc_event SET coinc_def_id = (SELECT new FROM _idmap_ WHERE old == coinc_def_id);
DELETE FROM coinc_definer WHERE coinc_def_id IN (SELECT old FROM _idmap_ WHERE old != new);

DROP INDEX tmpindex;
DROP TABLE _idmap_;

--
-- time_slide clean up
--

CREATE TEMPORARY TABLE _idmap_ AS
	SELECT
		time_slide_id AS old,
		(SELECT group_concat(instrument || "=" || offset) FROM time_slide AS time_slide_a WHERE time_slide_a.time_slide_id == time_slide.time_slide_id ORDER BY instrument) AS repr,
		NULL AS new
	FROM
		time_slide
	GROUP BY
		time_slide_id;
CREATE INDEX tmpindex ON _idmap_ (repr, old);

UPDATE _idmap_ SET new = (SELECT MIN(old) FROM _idmap_ AS a WHERE a.repr == _idmap_.repr);
DROP INDEX tmpindex;
CREATE INDEX tmpindex ON _idmap_ (old);

UPDATE coinc_event SET time_slide_id = (SELECT _idmap_.new FROM _idmap_ WHERE _idmap_.old == time_slide_id);
DELETE FROM time_slide WHERE time_slide_id IN (SELECT old FROM _idmap_ WHERE old != new);

DROP INDEX tmpindex;
DROP TABLE _idmap_;

-- VACUUM;
