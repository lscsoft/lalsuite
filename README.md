# LVK Algorithm Library Suite - LALSuite

LALSuite is a collection of data analysis routines whose primary application is
to search for and characterize astrophysical signals in gravitational-wave time
series data, in particular from ground-based detectors such as [LIGO], [Virgo],
and [KAGRA].

LALSuite is primarily written in C following the ISO/IEC 9899:1999 standard,
more commonly referred to as C99. It provides libraries written in C99, with
language bindings available for Python and Octave. It also provides command-line
applications written in C99 and Python.

LALSuite is licensed under the [GNU General Public License v2.0 or
later][GPL-2.0-or-later].

* [Project homepage][projectrepo]
* [Installation instructions][install]
* [Frequently asked questions][lalsuitefaq]
* [Reference documentation][lalsuitedocs]
* [Building from source][buildfromsrc]
* [Contributing guide][contributing]

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

If you make use of the [Python/Octave interfaces to LALSuite][swiglal], please
cite the following paper:

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

In addition, some codes contained in LALSuite may be directly based on one or
several scientific papers, which should be cited when using those specific
codes; some of these can be discovered through the documentation.

## Components

LALSuite is comprised of the following components:

- LAL: Core gravitational wave analysis routines
- LALFrame: LAL wrapping of the LIGO/Virgo Frame library
- LALMetaIO: LAL wrapping of the MetaIO LIGO_LW XML library
- LALSimulation: LAL routines for gravitational waveform and noise generation
- LALBurst: LAL routines for burst gravitational wave data analysis
- LALInspiral: LAL routines for inspiral and ringdown CBC gravitational wave
  data analysis
- LALInference: LAL routines for Bayesian inference data analysis
- LALPulsar: LAL routines for pulsar and continuous wave gravitational wave data
  analysis
- LALApps: Collection of gravitational wave data analysis codes and
  pipelines utilising the LAL libraries

## Project Librarians

- Adam Mercer
- Duncan Macleod
- Karl Wette

The librarians may be contacted via the [help desk][helpdesk].

## Documentation

[Installation instructions][install] and answers to [frequently asked questions][lalsuitefaq] may be found on the LALSuite [wiki][lalsuitewiki].

The [reference documentation][lalsuitedocs] for all LALSuite components is
generated using Doxygen. Documentation is available for the latest development
version of LALSuite, and for recent releases.

The *LAL Specification and Style Guide* document (commonly known as the *LAL
Spec*) contains the software specifications that code written for LALSuite
should conform to. The most recent version of the *LAL Spec* can be found
[here][lalspec].

## Releases

The latest release is LALSuite 7.26, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalsuite-7.26.tar.gz)), which is comprised of:

- LAL 7.7.0, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lal-7.7.0.tar.xz))
- LALFrame 3.0.7, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalframe-3.0.7.tar.xz))
- LALMetaIO 4.0.6, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalmetaio-4.0.6.tar.xz))
- LALSimulation 6.2.0, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalsimulation-6.2.0.tar.xz))
- LALBurst 2.0.7, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalburst-2.0.7.tar.xz))
- LALInspiral 5.0.3, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalinspiral-5.0.3.tar.xz))
- LALInference 4.1.9, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalinference-4.1.9.tar.xz))
- LALPulsar 7.1.1, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalpulsar-7.1.1.tar.xz))
- LALApps 10.1.0, released 22 May 2025 ([source tarball](https://software.igwn.org/lscsoft/source/lalsuite/lalapps-10.1.0.tar.xz))

Previous releases can be found by visiting the [IGWN analysis software source
code server][lalsuitesources].

## Sources

The LALSuite sources are located in the [LALSuite Git repository][projectrepo].
See the following [clonerepo] on how to get started.

## Mailing Lists

There are several mailing lists related to LALSuite development:

- *LAL Announce*: General LALSuite announcements, including releases
    ([subscribe/unsubscribe][lal-announce-s], [archive][lal-announce-a])
- *LAL Discuss*: Discussions of modifying/using LALSuite
    ([subscribe/unsubscribe][lal-discuss-s], [archive][lal-discuss-a])
- *Computing Discuss*: General IGWN Computing discussion
    ([subscribe/unsubscribe][cmp-discuss-s], [archive][cmp-discuss-a])

Developers should subscribe to the *LAL Discuss* and *Computing Discuss* mailing
lists as important development information is regularly posted and discussed on
these lists.

## Reporting Issues

If you have `ligo.org` authentication, please report issues directly through
GitLab. Otherwise, you can use the [help desk][helpdesk] to send bug reports by
e-mail.

Before you file a ticket, please read and search through the list of
[current][currentissues] and [previous][previousissues] issues to determine if
your bug has already been reported.

If an issue already exists and has not been fixed, add any additional
information to the existing report. If your bug exists and has been fixed,
upgrade to the version detailed in the issue to confirm if it has been fixed
correctly. If it was not, please reopen the issue.

Please include as much detail as possible to reproduce the
error, including information about your operating system and the version of each
(relevant) component of LALSuite. If possible, please include a brief,
self-contained code example that demonstrates the problem.

Note that when an issue is marked as *Confidential*, currently this means that
most internal users will also not be able to see it, but only a small number of
people with reporter, developer or maintainer status.

## Contributing

The [guide to contributing to LALSuite][contributing] explains how to contribute
fixes or new features using the fork and pull workflow. Please read and follow
these directions.

[GPL-2.0-or-later]: https://spdx.org/licenses/GPL-2.0-or-later.html
[KAGRA]:            https://gwcenter.icrr.u-tokyo.ac.jp/en/
[LIGO]:             https://www.ligo.org
[Virgo]:            https://www.virgo-gw.eu/
[buildfromsrc]:     https://git.ligo.org/lscsoft/lalsuite/-/wikis/BUILD
[clonerepo]:        https://git.ligo.org/lscsoft/lalsuite/-/wikis/BUILD#cloning-the-git-repository
[cmp-discuss-a]:    https://sympa.ligo.org/wws/info/computing-discuss
[cmp-discuss-s]:    https://grouper.ligo.org/selfmanage/computing-discuss
[contributing]:     https://git.ligo.org/lscsoft/lalsuite/-/blob/master/CONTRIBUTING.md
[currentissues]:    https://git.ligo.org/lscsoft/lalsuite/issues
[doi]:              https://doi.org/10.7935/GT1W-FZ16
[helpdesk]:         mailto:contact+lscsoft-lalsuite-1438-issue-@support.ligo.org
[install]:          https://git.ligo.org/lscsoft/lalsuite/-/wikis/INSTALL
[lal-announce-a]:   https://sympa.ligo.org/wws/info/lal-announce
[lal-announce-s]:   https://grouper.ligo.org/selfmanage/lal-announce
[lal-discuss-a]:    https://sympa.ligo.org/wws/info/lal-discuss
[lal-discuss-s]:    https://grouper.ligo.org/selfmanage/lal-discuss
[lalspec]:          https://dcc.ligo.org/LIGO-T990030/public
[lalsuitedocs]:     https://docs.ligo.org/lscsoft/lalsuite
[lalsuitefaq]:      https://git.ligo.org/lscsoft/lalsuite/-/wikis/FAQ
[lalsuitesources]:  http://software.ligo.org/lscsoft/source/lalsuite/?C=M;O=D
[lalsuitewiki]:     https://git.ligo.org/lscsoft/lalsuite/-/wikis
[previousissues]:   https://git.ligo.org/lscsoft/lalsuite/issues?scope=all&utf8=%E2%9C%93&state=closed
[projectrepo]:      https://git.ligo.org/lscsoft/lalsuite
[swiglal]:          https://docs.ligo.org/lscsoft/lalsuite/dev/swiglal_tutorial.html
