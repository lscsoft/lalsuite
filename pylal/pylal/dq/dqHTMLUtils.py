#!/usr/bin/env python

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import re,socket
from glue import markup

from glue import git_version

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides a few extensions to glue.markup to streamline GEO/LIGO detector characterisation tools that use very similar web interfaces
"""

# =============================================================================
# Write table
# =============================================================================

def write_table(page, headers, data, cl=''):

  """
    Write table into glue.markup.page object. headers are written with <th>,
    multiple columns of data are written with <td>. Classes include:
      * "", default class writes standard table with single column of headers
            and multiple columns of data
      * "list", writes table of 'header[i]: data[i]' definition-style entries

    Arguments:

      page : glue.markup.page
        page object into which to write table
      headers : list
        list of table header elements
      data : list
        list (or nested list) of table data elements, list of lists used for
        multiple rows
 
    Keyword arguments:

      cl : string
        name for HTML table class, cl='list' treats special case above

  """

  # open table
  page.table(class_=cl)

  # list: print two vertical columns of header:data pairs
  if cl=='list':
    for i in range(len(headers)):

      page.tr()
      page.th(str(headers[i]))
      page.td(str(data[i]))
      page.tr.close()

  # otherwise print 'standard' table with single header row and multiple data
  # rows
  else:
    page.tr()
    if len(headers)==1:
      page.th(str(headers[0]), colspan="100%")
    else:
      for n in headers:
        page.th(str(n))
    page.tr.close()

    if data and not re.search('list',str(type(data[0]))):
      data = [data]

    for row in data:
      page.tr()
      for item in map(str, row):
        page.td(item)

  page.table.close()

  return page

# =============================================================================
# Write <head>
# =============================================================================

def write_head(title, css, js, base=None, refresh=None, jquery=True):

  """
    Returns glue.markup.page object with <head> tag filled.

    Arguments:

      title : string
        text for <title> tag
      css : string
        relative path to style sheet
      js : string
        relative path to javascript

    Keyword arguments:

      base : string
        absolute http(s) path of url base
      refresh : int
        number of seconds after which to refresh page automatically
      jquery : [ True | False ]
        import jquery AJAX script in header, default: True
  """

  # generate object
  page = markup.page(mode="strict_html")
  page._escape = False

  # open head
  page.head()
  # add base
  if base:
    page.base(href=base)
  # add html auto-refresh
  if refresh:
    page.meta(http_equiv="refresh", content="%s" % refresh)
  # link stylesheet
  if isinstance(css, str): css = [css]
  for c in css:
    page.link(media="all", href=c, type="text/css", rel="stylesheet")
  # add title
  page.title(title)

  if jquery:
    page.script("", src="http://ajax.googleapis.com/ajax/libs/jquery/1.7.0"\
                    "/jquery.min.js", type="text/javascript")
  if isinstance(js, str): js = [js]
  for j in js:
    page.script("", src=j, type="text/javascript")
  page.head.close()

  return page

# =============================================================================
# Write <div id="header">
# =============================================================================

def write_banner(title, text=""):

  """
    Returns glue.markup.page object for <div id="header">
  """

  page = markup.page(mode="strict_html")
  page._escape = False

  page.div(class_="content", id="header")
  page.div()
  page.h1(title)
  page.h3(text)
  page.div.close()

  page.div.close()

  return page

# =============================================================================
# Write <div id="menubar">
# =============================================================================

def write_menu(sections, pages, current=None):

  """
    Returns glue.markup.page object for <div id="menubar">, constructing menu
    in HTML.

    Arguments:

      sections : list
        ordered list of menu entry names
      pages : dict
        dict of section:href pairs holding link paths for each element of
        sections list
  """

  page = markup.page(mode="strict_html")
  page._escape = False

  page.div(id="menubar")

  for i,sec in enumerate(sections):
    if sec==current: cl = "menulink selected"
    else:            cl = "menulink"
    page.a(sec, id="link_%d" % i, class_=cl, href=pages[sec])

  page.script("",type="text/javascript")

  page.div.close()

  return page

# =============================================================================
# Initialise page
# =============================================================================

def init_page(head, banner, menu, **kwargs):

  """
    Initialise html into markup page, including <head> tag, banner and menu.
    Pass further html elements to the body tag using the kwargs.
  """

  # write html
  page = markup.page()
  page._escape = False

  # initialise page
  page.add("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" "+\
           "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">")
  page.html(xmlns="http://www.w3.org/1999/xhtml", lang="en",\
            **{"xml:lang":"en"})
  page.add(head())

  # open body
  page.body(**kwargs)

  # open container for page (needed to position footer)
  page.div(id="container")
  # add banner
  page.add(banner())
  # open content (tab below banner and above footer)
  page.div(id="content")
  page.div()
  # print menu
  page.add(menu())
  # initialise maintab
  page.div(id="maintab")

  return page

# =============================================================================
# Close page
# =============================================================================

def close_page(page, footer=False):

  """
    Close content, maintab, and container divs, write footer and close
    <body> and <html> tags.
  """

  # close maintab
  page.div.close()
  # close content tab
  page.div.close()
  page.div.close()
  # close container
  page.div.close()
  # add footer
  if footer:
    page.add(footer())

  # close page
  page.body.close()
  page.html.close()

  return page

# =============================================================================
# Write glossary
# =============================================================================

def write_glossary(page, terms):

  """
    Write a glossary of DQ terms into the glue.markup.page object page using the
    list of (term,definition) tuples terms.
  """

  page.h2()
  page.add('Glossary')
  page.h2.close()
  page.add('This section gives a glossary of terms used in this analysis')

  i=1

  # write section
  page = write_h(page, 'Glossary', [i], cl=3)
  page.div(id="div_%s" % (i), style="display: none;")

  page.add('The LIGO-Virgo acronym wiki can be found on')
  page.a(href="https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/Acronyms")
  page.add('this page')
  page.a.close()

  # write glossary table
  tmp = {}
  for (term, text) in terms:
    tmp[term.title()] = text
  th = sorted(tmp.keys())
  td = [tmp[t] for t in th]
  page = write_table(page, th, td, 'list')

  page.div.close()

  return page

# =============================================================================
# Write <hX>
# =============================================================================

def write_h(page, title, id, cl=3, toggle=True):

  """
    Write hX header into glue.markup.page object page with toggling. Text
    contained within title is printed while link is constructed and toggled
    using javascipt toggleVisible function
  """

  if not (isinstance(id, list) or isinstance(id, tuple)):
    id = [id]
  id = list(map(str, id))

  kwargs = {}
  if toggle:
    kwargs['onclick'] = "toggleVisible('%s');" % '.'.join(id)

  page.input(id="input_%s" % '.'.join(id), type="button", class_="h%s" % cl,\
             value=title, **kwargs)

  return page

# =============================================================================
# Link image
# =============================================================================

def link_image(page, href, src, alt, title, **kwargs):

  """
    Link image into glue.markup.page object page with standard options
  """

  page.a(href=href, title=title, **kwargs)
  page.img(src=src, alt=alt, **kwargs)
  page.a.close()

  return page

# =============================================================================
# Link file
# =============================================================================

def link_file(page, href, text):

  """
    Link file into glue.markup.page object page, with associated text and
    options.
  """

  page.a(text, href=href, rel="external")

  return page

# =============================================================================
# Construct HTML base
# =============================================================================

def get_ldas_url():
  """
    Returns the url for http access to this host, based on its domain name
    Returns None-type if you're not on a recognised LIGO-Virgo cluster.
    Should not be used when not on an LDAS cluster, can't tell the difference 
    between nodes on CDS network at a site, and those on LDAS, for example...

    Example:
    >>> socket.getfqdn()
    ldas-pcdev1.ligo-la.caltech.edu
    >>> dqHTMLUtils.get_ldas_url()
    'https://ldas-jobs.ligo-la.caltech.edu'
  """

  # get fully qualified domain name (FQDN)
  fqdn = socket.getfqdn()

  # work out web root
  if re.search('ligo.caltech', fqdn):
    root = "https://ldas-jobs.ligo.caltech.edu"
  elif re.search('ligo-wa', fqdn):
    root = "https://ldas-jobs.ligo-wa.caltech.edu"
  elif re.search('ligo-la', fqdn):
    root = "https://ldas-jobs.ligo-la.caltech.edu"
  elif re.search('atlas', fqdn):
    match = re.search('atlas[13]', fqdn)
    if match:
      host = match.group()
      root = "https://%s.atlas.aei.uni-hannover.de" % host
    else:
      match = "https://atlas1.atlas.aei.uni-hannover.de"
  elif re.search('phy.syr', fqdn):
    root = "https://sugar-jobs.phy.syr.edu"
  elif re.search('astro.cf', fqdn):
    root = "https://gravity.astro.cf.ac.uk"
  elif re.search('phys.uwm', fqdn):
    root = "https://ldas-jobs.phys.uwm.edu"
  else:
    root = None

  return root
