# Short Fourier Transforms data files

This directory contains Short Fourier Transform (SFT) data files. SFTs are a
standard input data source for continuous gravitational wave analysis
pipelines. An SFT is stored in an SFT file in a binary format. Each SFT file may
contain multiple SFTs from the same detector. The SFT file format and naming
scheme are further described in the [SFT specification][SFTspec].

The SFTs in this directory are organised as follows.

## Common prefix

SFTs are organised into directories with a common prefix, and a
broadband/narrowband-specific suffix. The common prefix is named according to
the following scheme:

     {common-prefix} = {det}_{Tsft}SFT
                       _O{run}{kind}+R{rev}+C{chan}+W{win}
                       _{bandw}

(Note that the right hand side is broken into separate lines for readability; no
whitespace appears in the actual directory names.)

Terms in braces denote the following components:

- `det`: Gravitational wave detector of the data. Common choices are:
  - `H1`: LIGO Hanford (4 km)
  - `K1`: KAGRA
  - `L1`: LIGO Livingston
  - `V1`: Virgo

- `Tsft`: Timespan of the data, in seconds, which is Fourier transformed to
          create each SFT. A typical choice is 1800 seconds, which is generally
          appropriate for searches for continuous gravitational waves from
          isolated neutron stars.

- `run`: Observing run number. For example, `run = 1` corresponds to the first
         observing run (O1) of the LIGO Observatory, from September 12, 2015
         through January 19, 2016.

- `rev`: Revision number of the SFT dataset. Typically `rev = 1` the first time
         the dataset is generated. If, however, it is later found that a dataset
         was generated incorrectly, it may be regenerated with an incremented
         revision number, e.g. `rev = 2`. Please check carefully that you are
         using the latest revision of an SFT dataset!

- `chan`: Name of the data channel (without the detector prefix) from which the
          SFTs are generated. Due to limitations on the SFT file/directory
          naming scheme, the channel name is stripped of any non-alphanumeric
          characters, which may make it difficult to read. For example, `chan =
          GDSCALIBSTRAINCLEANGATEDG01` corresponds to the data channel
          `{det}:GDS-CALIB_STRAIN_CLEAN_GATED_G01`.

- `win`: Windowing applied to the SFTs. Typical choices are:

  - `TKEY{n}`: Tukey windowing, with parameter `beta = n / 5000`. A standard
               choice for `beta` is 0.001, which corresponds to `n = 5`. Note
               that, with this choice of `beta`, the Tukey window is very close
               to a rectangular window. This windowing is chosen to mitigate any
               transients at the start and end of the SFT. The loss in signal
               power due to the windowing is on the order of 0.1%; see [footnote
               65 in this paper][S5CasA].

  - `HANN`: Hann windowing.

- `bandw`: Bandwidth of the SFTs. This may be one of two choices:

  - `BROADBAND`: SFTs which are broadband in frequency, typically covering
                 the entire useful/sensitive band of the detectors.
                 Conversely, these SFT files are typically limited in time;
                 each SFT file contains only 1 SFT, and spans a time `Tsft`.

  - `NARROWBAND`: SFTs which are narrowband in frequency, typically covering
                  only a few Hz. Conversely, these SFT files typically cover
                  longer time ranges; each SFT file may contain more than one
                  SFT, and each SFT file may therefore span a time much
                  greater than `Tsft`.

## Broadband SFT suffix

Broadband SFT files (where `bandw = BROADBAND`) are organised into directories
based on their starting GPS time, as follows:

     {broadband-sft-path} = {common-prefix}-{start-1e6}/
                            {broadband-sft-filename}

where `start-1e6` is the starting GPS time of the file divided by 1e6, rounded
down. For example, an SFT starting at GPS time 1,389,258,149 would be found in
the `-1389/` directory.

Within each directory, SFT files are named as follows:

     {broadband-sft-filename} = {site}-{num}
                                _{det}_{Tsft}SFT
                                _O{run}{kind}+R{rev}+C{chan}+W{win}
                                -{start}-{span}.sft

where, in addition to the top-level directory components:

- `site`: Gravitational wave observatory site. This is always the first letter
          of `det`, e.g. `H`, `K`, `L`, `V`.
- `num`: Number of SFTs contained in the file.
- `start`: Starting GPS time of the SFT data in the file.
- `span`: Total timespan of all the SFTs in the file. If `num = 1`, then
          `span = Tsft`, otherwise `span > Tsft`.

## Narrowband SFT subdirectories

Narrowband SFT files (where `bandw = NARROWBAND`) are organised into directories
based on their frequency bands, as follows:

     {narrowband-sft-path} = {common-prefix}-NBF{freq-Hz}Hz{freq-bin}W{width-Hz}Hz{width-bin}/
                             {narrowband-sft-filename}

`{freq-Hz}` and `{freq-Hz}` are the result of the integer division of the
starting bin of the SFT by `Tsft`, and the remainder of that division,
respectively. For example, a `Tsft = 1800` narrowband SFT starting at bin
180,900 would have

     freq-Hz = 180,900 / 1800 s = 100.5 Hz = 100 Hz, rounded down
     width-Hz = 180,900 - 100 Hz * 1800 s = 900 bins
     i.e.
     starting bin = 100 Hz * 1800 s + 900 bins = bin 180,900

Similarly, `{width-Hz}` and `{width-Hz}` are calculated from the bandwidth of
the SFT.

Within each directory, SFT files are named as follows:

     narrowband-sft-filename = {site}-{num}
                               _{det}_{Tsft}SFT
                               _O{run}{kind}+R{rev}+C{chan}+W{win}
                               _NBF{freq-Hz}Hz{freq-bin}W{width-Hz}Hz{width-bin}
                               -{start}-{span}.sft

where the components are as described above.

## For further information

Please consult the [SFT specification][SFTspec].

[SFTspec]: https://dcc.ligo.org/LIGO-T040164-v4/public
[S5CasA]:  http://doi.org/10.1088/0004-637X/722/2/1504
