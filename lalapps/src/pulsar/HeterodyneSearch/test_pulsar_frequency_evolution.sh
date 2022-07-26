# generate frequency evolution of pulsar given by test.par

grep . test.par

lalapps_pulsar_frequency_evolution --detector L1 --par-file test.par --start 1100000000 --timespan 3600 --deltat 60 --output-dir .

# compare to reference result

paste ./frequency_evolution_L1.txt ./reference/frequency_evolution_L1.txt | awk '
function relative_error(col, x, x0) {
   relerr = sqrt( (x - x0)^2 / x0^2 )
   if (relerr > 1e-6) {
      printf "line %i, column %s : %g does not compare to reference %g\n", FNR, col, x, x0
      exit(1)
   }
}
{
   relative_error("1 (gps)",       $1, $6)
   relative_error("2 (freq)",      $2, $7)
   relative_error("3 (dopplerss)", $3, $8)
   relative_error("4 (dopplerbs)", $4, $9)
   relative_error("5 (doppler)",   $5, $10)
}
'
