# LALInspiral

LAL routines for inspiral and ringdown CBC gravitational wave data
analysis.

## Acknowledgment

We request that any academic report, publication, or other academic
disclosure of results derived from the use of this software acknowledge
the use of the software by an appropriate acknowledgment or citation.

The LALSuite software suite as a whole may be linked to using the DOI
[10.7935/GT1W-FZ16][doi], or cited using the following BibTeX entry:

    @misc{lalsuite,
        author         = "{LIGO Scientific Collaboration}
                          and {Virgo Collaboration}
                          and {KAGRA Collaboration}",
        title          = "{LVK} {A}lgorithm {L}ibrary - {LALS}uite",
        howpublished   = "Free software (GPL)",
        doi            = "10.7935/GT1W-FZ16",
        year           = "2018"
    }

If you make use of the [Python/Octave interfaces to LALInspiral][swiglal],
please cite the following paper:

    @article{swiglal,
        title     = "{SWIGLAL}: {P}ython and {O}ctave interfaces to the
                     {LALSuite} gravitational-wave data analysis libraries}",
        author    = "Karl Wette",
        journal   = "SoftwareX",
        volume    = "12",
        pages     = "100634",
        year      = "2020",
        doi       = "10.1016/j.softx.2020.100634"
    }

In addition, some codes contained in LALInspiral may be directly based on one or
several scientific papers, which should be cited when using those specific
codes; some of these can be discovered through the documentation.

## Basic Build Instructions

     ./configure --prefix=...
     make
     make install

Please read the [LALSuite install how-to][install] for more detailed
build instructions.

## Bug Reporting

Please visit the [LALSuite bug reporting system][bugs] (LIGO.org
authentication required to submit new issues), or use the [e-mail
helpdesk][helpdesk].

## For More Information

See the [LALSuite README][readme].

[bugs]:             https://git.ligo.org/lscsoft/lalsuite/issues/
[doi]:              https://doi.org/10.7935/GT1W-FZ16
[helpdesk]:         mailto:contact+lscsoft-lalsuite-1438-issue-@support.ligo.org
[install]:          https://git.ligo.org/lscsoft/lalsuite/-/wikis/INSTALL
[readme]:           https://git.ligo.org/lscsoft/lalsuite/-/blob/master/README.md
[swiglal]:          https://lscsoft.docs.ligo.org/lalsuite/dev/swiglal_tutorial.html
