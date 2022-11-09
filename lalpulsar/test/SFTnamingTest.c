/*
 * Copyright (C) 2022 Karl Wette
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

/**
 * \file
 * \ingroup SFTfileIO_h
 * \author K. Wette
 *
 * \brief Test the SFT file naming routines
 */

#include <lal/SFTfileIO.h>

int main(void) {

  const char* SFT_filenames[] = {
    // A 30-minute private H2 SFT
    "/path/to/H-1_H2_1800SFT-735627918-1800.sft",
    // A 30-minute private L1 SFT with the description "S2"
    "L-1_L1_1800SFT_S2-733467931-1800.sft",
    // A private file with 1887 30-minute H1 SFTs with gaps in time
    "H-1887_H1_1800SFT-733467931-4622400.sft",
    // A 60-second private GEO SFT with the description "S3hot"
    "G-1_G1_60SFT_S3hot-732465218-60.sft",
    // 5 public broad-band H1 SFTs from the O2 observing run; revision 1 of the public SFTs were
    // created from the H1:DCH-CLEAN_STRAIN_C02 channel, Tukey-windowed with parameter = 0.001
    "H-1_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5-1257800000-1800.sft",
    "H-1_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5-1257801800-1800.sft",
    "H-1_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5-1257803600-1800.sft",
    "H-1_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5-1257805400-1800.sft",
    "H-1_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5-1257807200-1800.sft",
    // The equivalent narrow-band SFTs over the frequency range 10–95.5 Hz, each SFT
    // containing 8 Hz of data except the last SFT which contains 5.5 Hz of data
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0010Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0018Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0026Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0034Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0042Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0050Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0058Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0066Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0074Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0082Hz0W0008Hz0-1257800000-9000.sft",
    "H-5_H1_1800SFT_O2RUN+R1+CDCHCLEANSTRAINC02+WTKEY5_NBF0090Hz0W0005Hz900-1257800000-9000.sft",
  };

  const SFTFilenameSpec SFT_spec[] = {
    // A 30-minute private H2 SFT
    { .path = "/path/to", .numSFTs = 1, .detector = "H2", .SFTtimebase = 1800, .gpsStart = 735627918, .SFTspan = 1800 },
    // A 30-minute private L1 SFT with the description "S2"
    { .numSFTs = 1, .detector = "L1", .SFTtimebase = 1800, .privMisc = "S2", .gpsStart = 733467931, .SFTspan = 1800 },
    // A private file with 1887 30-minute H1 SFTs with gaps in time
    { .numSFTs = 1887, .detector = "H1", .SFTtimebase = 1800, .gpsStart = 733467931, .SFTspan = 4622400 },
    // A 60-second private GEO SFT with the description "S3hot"
    { .numSFTs = 1, .detector = "G1", .SFTtimebase = 60, .privMisc = "S3hot", .gpsStart = 732465218, .SFTspan = 60 },
    // 5 public broad-band H1 SFTs from the O2 observing run; revision 1 of the public SFTs were created from the H1:DCH-CLEAN_STRAIN_C02 channel, Tukey-windowed with parameter = 0.001
    { .numSFTs = 1, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 1800,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02" },
    { .numSFTs = 1, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257801800, .SFTspan = 1800,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02" },
    { .numSFTs = 1, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257803600, .SFTspan = 1800,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02" },
    { .numSFTs = 1, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257805400, .SFTspan = 1800,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02" },
    { .numSFTs = 1, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257807200, .SFTspan = 1800,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02" },
    // The equivalent narrow-band SFTs over the frequency range 10–95 Hz, each SFT containing 8 Hz of data except the last SFT which contains 5 Hz of data
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 10, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 18, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 26, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 34, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 42, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 50, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 58, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 66, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 74, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 82, .nbFirstBinRem = 0, .nbBinWidthFreq = 8, .nbBinWidthRem = 0 },
    { .numSFTs = 5, .detector = "H1", .SFTtimebase = 1800, .window_type = "tukey", .window_param = 0.001, .gpsStart = 1257800000, .SFTspan = 9000,
      .pubObsRun = 2, .pubObsKind = "RUN", .pubRevision = 1, .pubChannel = "H1:DCH-CLEAN_STRAIN_C02",
      .nbFirstBinFreq = 90, .nbFirstBinRem = 0, .nbBinWidthFreq = 5, .nbBinWidthRem = 900 },
  };

  for (size_t i = 0; i < XLAL_NUM_ELEM(SFT_spec); ++i) {
    char *fname = XLALBuildSFTFilenameFromSpec(&SFT_spec[i]);
    XLAL_CHECK_MAIN( fname != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN( strcmp( fname, SFT_filenames[i] ) == 0, XLAL_EFAILED, "SFT filename '%s' should be '%s'", fname, SFT_filenames[i] );
    XLALFree( fname );
  }

  for (size_t i = 0; i < XLAL_NUM_ELEM(SFT_spec); ++i) {
    SFTFilenameSpec spec;
    XLAL_CHECK_MAIN( XLALParseSFTFilenameIntoSpec(&spec, SFT_filenames[i]) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( strcmp( spec.path, SFT_spec[i].path ) == 0, XLAL_EFAILED, "path '%s' should be '%s'", spec.path, SFT_spec[i].path );
    XLAL_CHECK_MAIN( spec.numSFTs == SFT_spec[i].numSFTs, XLAL_EFAILED, "numSFTs '%d' should be '%d'", spec.numSFTs, SFT_spec[i].numSFTs );
    XLAL_CHECK_MAIN( strcmp( spec.detector, SFT_spec[i].detector ) == 0, XLAL_EFAILED, "detector '%s' should be '%s'", spec.detector, SFT_spec[i].detector );
    XLAL_CHECK_MAIN( spec.SFTtimebase == SFT_spec[i].SFTtimebase, XLAL_EFAILED, "SFTtimebase '%d' should be '%d'", spec.SFTtimebase, SFT_spec[i].SFTtimebase );
    XLAL_CHECK_MAIN( strcmp( spec.window_type, SFT_spec[i].window_type ) == 0, XLAL_EFAILED, "window_type '%s' should be '%s'", spec.window_type, SFT_spec[i].window_type );
    XLAL_CHECK_MAIN( spec.window_param == SFT_spec[i].window_param, XLAL_EFAILED, "window_param '%g' should be '%g'", spec.window_param, SFT_spec[i].window_param );
    XLAL_CHECK_MAIN( spec.gpsStart == SFT_spec[i].gpsStart, XLAL_EFAILED, "gpsStart '%d' should be '%d'", spec.gpsStart, SFT_spec[i].gpsStart );
    XLAL_CHECK_MAIN( spec.SFTspan == SFT_spec[i].SFTspan, XLAL_EFAILED, "SFTspan '%d' should be '%d'", spec.SFTspan, SFT_spec[i].SFTspan );
    XLAL_CHECK_MAIN( strcmp( spec.privMisc, SFT_spec[i].privMisc ) == 0, XLAL_EFAILED, "privMisc '%s' should be '%s'", spec.privMisc, SFT_spec[i].privMisc );
    XLAL_CHECK_MAIN( spec.pubObsRun == SFT_spec[i].pubObsRun, XLAL_EFAILED, "pubObsRun '%d' should be '%d'", spec.pubObsRun, SFT_spec[i].pubObsRun );
    XLAL_CHECK_MAIN( strcmp( spec.pubObsKind, SFT_spec[i].pubObsKind ) == 0, XLAL_EFAILED, "pubObsKind '%s' should be '%s'", spec.pubObsKind, SFT_spec[i].pubObsKind );
    XLAL_CHECK_MAIN( spec.pubRevision == SFT_spec[i].pubRevision, XLAL_EFAILED, "pubRevision '%d' should be '%d'", spec.pubRevision, SFT_spec[i].pubRevision );
    if (strlen(spec.pubChannel) > 0) {
      const char *SFT_spec_pubChannel_stripped = "DCHCLEANSTRAINC02";
      XLAL_CHECK_MAIN( strcmp( spec.pubChannel, SFT_spec_pubChannel_stripped ) == 0, XLAL_EFAILED, "pubChannel '%s' should be '%s'", spec.pubChannel, SFT_spec_pubChannel_stripped );
    }
    XLAL_CHECK_MAIN( spec.nbFirstBinFreq == SFT_spec[i].nbFirstBinFreq, XLAL_EFAILED, "nbFirstBinFreq '%d' should be '%d'", spec.nbFirstBinFreq, SFT_spec[i].nbFirstBinFreq );
    XLAL_CHECK_MAIN( spec.nbFirstBinRem == SFT_spec[i].nbFirstBinRem, XLAL_EFAILED, "nbFirstBinRem '%d' should be '%d'", spec.nbFirstBinRem, SFT_spec[i].nbFirstBinRem );
    XLAL_CHECK_MAIN( spec.nbBinWidthFreq == SFT_spec[i].nbBinWidthFreq, XLAL_EFAILED, "nbBinWidthFreq '%d' should be '%d'", spec.nbBinWidthFreq, SFT_spec[i].nbBinWidthFreq );
    XLAL_CHECK_MAIN( spec.nbBinWidthRem == SFT_spec[i].nbBinWidthRem, XLAL_EFAILED, "nbBinWidthRem '%d' should be '%d'", spec.nbBinWidthRem, SFT_spec[i].nbBinWidthRem );
  }

  LALCheckMemoryLeaks();

}
