CREATE TEMPORARY TABLE coinc_inspiral_max AS SELECT max(snr) AS maxsnr, end_time AS end_time FROM coinc_inspiral GROUP BY end_time;

SELECT * FROM coinc_inspiral JOIN coinc_inspiral_max ON (coinc_inspiral.snr == coinc_inspiral_max.maxsnr and coinc_inspiral_max.end_time == coinc_inspiral.end_time);
