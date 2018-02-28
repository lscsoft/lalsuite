"""
lightweight database dumper
Copyright (C) 2003 Duncan Brown
This file is part of the lightweight datapase dumper (ldbd)

The ldbd module provides classes for manipulating LIGO metadata database
tables.

References:
http://www.ligo.caltech.edu/docs/T/T990101-02.pdf
http://www.ligo.caltech.edu/docs/T/T990023-01.pdf
http://ldas-sw.ligo.caltech.edu/doc/db2/doc/html/ilwdformat.html

This file is part of the Grid LSC User Environment (GLUE)

GLUE is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from glue import git_version
__author__ = 'Duncan Brown <dbrown@ligo.caltech.edu>'
__date__ = git_version.date 
__version__ = git_version.id

import os
import sys
import string
import re
import csv
import exceptions
try:
  import DB2
except:
  pass

from glue.ligolw.types import string_format_func

"""
create the csv parser and initialize a dialect for LIGO_LW streams
"""

class LIGOLWStream(csv.Dialect):
  """
  Create a csv parser dialect for parsing LIGO_LW streams
  """
  delimiter = ','
  doublequote = False
  escapechar = '\\'
  lineterminator = '\n'
  quotechar = '"'
  quoting = csv.QUOTE_ALL
  skipinitialspace = True

csv.register_dialect("LIGOLWStream",LIGOLWStream)


class LIGOLwParseError(exceptions.Exception):
  """Error parsing LIGO lightweight XML file"""
  pass


class LIGOLwDBError(exceptions.Exception):
  """Error interacting with database"""
  pass


class Xlator(dict):
  """
  All in one multiple string substitution class from the python cookbook
  """
  def _make_regex(self):
    """
    Build a re object based on keys in the current dictionary
    """
    return re.compile("|".join(map(re.escape, self.keys())))

  def __call__(self, match):
    """
    Handler invoked for each regex match
    """
    return self[match.group(0)]

  def xlat(self, text):
    """
    Translate text, returns the modified text
    """
    return self._make_regex().sub(self,text)


class LIGOMetadataDatabase:
  """
  Contains a tuple of tables in the order that insertions should
  ocour and a dictionary of mandatory uniquw id fields for each
  table that must be populated.
  """
  def __init__(self,database):
    """
    database = the name of the LIGO database to initalize
    """
    self.database = database
    self.uniqueids = {}
    conn = DB2.connect(dsn=database, uid='', pwd='')
    curs = conn.cursor()
    curs.execute("SELECT tabname FROM syscat.tables WHERE definer<>'SYSIBM' "
      "AND TYPE='T' ORDER BY PARENTS ASC")
    self.tables = curs.fetchall()
    curs.execute("SELECT tabname, colname FROM syscat.columns " +
      "WHERE typename = 'CHARACTER' AND length = 13 AND nulls = 'N'")
    for tab, col in curs.fetchall():
      tab = tab.lower()
      col = col.lower()
      try:
        self.uniqueids[tab][col] = 'ilwd:char'
      except KeyError:
        self.uniqueids[tab] = {}
        self.uniqueids[tab][col] = 'ilwd:char'
    curs.close()
    conn.close()


class UniqueIds:
  """
  Contains a dictionary of unique ids which can be queried based
  on name. If a unique id does not exist in the dictionaty, one
  is fetched from the database. 
  """
  def __init__(self,curs):
    """
    curs = database cursor to the currently open database
    """
    self.uqids = {}
    self.curs = curs

  def lookup(self,istring):
    """
    istring = the ilwd:char string corresponding to a unique id
    """
    try:
      return self.uqids[istring]
    except KeyError:
      curs = self.curs
      curs.execute('VALUES BLOB(GENERATE_UNIQUE())')
      self.uqids[istring] = curs.fetchone()[0]
      return self.uqids[istring]


class LIGOLwParser:
  """
  Provides methods for parsing the data from a LIGO lightweight XML
  file parsed with pyRXP into a dictionary
  """

  def __init__(self):
    """
    Initializes a LIGO lightweight XML parser with the necessary 
    regular expressions and function for tuple translation
    """
    self.tabrx = re.compile(r'(\A[a-z0-9_]+:|\A)([a-z0-9_]+):table\Z')
    self.colrx = re.compile(r'(\A[a-z0-9_]+:|\A)([a-z0-9_]+:|\A)([a-z0-9_]+)\Z')
    self.llsrx = re.compile(r'\A\s*"')
    self.rlsrx = re.compile(r'"\s*\Z')
    self.licrx = re.compile(r'\A\s+"')
    self.ricrx = re.compile(r'"*\s*\Z')
    self.octrx = re.compile(r'\A\\[0-9][0-9][0-9]')
    self.dlmrx = re.compile(r'\\,')
    self.unique = None
    self.types = {
      'int_2s' : int,
      'int_4s' : int,
      'int_8s' : long,
      'real_4' : float,
      'real_8' : float,
      'lstring' : self.__lstring,
      'ilwd:char' : self.__ilwdchar,
      'ilwd:char_u' : self.__ilwdchar
    }
    self.xmltostr = Xlator({ r'&amp;' : r'&', r'&gt;' : r'>', r'&lt;' : r'<','\\\\' : '\\'}) # Note: see https://www.gravity.phy.syr.edu/dokuwiki/doku.php?id=rpfisher:gluebughunt if this is confusing, the parser just cleanly handles the conversion of everything

  def __del__(self):
    if self.unique:
      del self.unique

  def __lstring(self,lstr):
    """
    Returns a parsed lstring by stripping out and instances of
    the escaped delimiter. Sometimes the raw lstring has whitespace
    and a double quote at the beginning or end. If present, these
    are removed.
    """
    lstr = self.llsrx.sub('',lstr.encode('ascii'))
    lstr = self.rlsrx.sub('',lstr)
    lstr = self.xmltostr.xlat(lstr)
    lstr = self.dlmrx.sub(',',lstr)
    return lstr
    
  def __ilwdchar(self,istr):
    """
    If the ilwd:char field contains octal data, it is translated
    to a binary string and returned. Otherwise a lookup is done
    in the unique id dictionary and a binary string containing the
    correct unique id is returned.
    """
    istr_orig = istr
    istr = self.licrx.sub('',istr.encode('ascii'))
    istr = self.ricrx.sub('',istr)
    if self.octrx.match(istr):
      exec "istr = '"+istr+"'"
      # if the DB2 module is loaded, the string should be converted
      # to an instance of the DB2.Binary class. If not, leave it as 
      # a string containing binary data.
      try:
        istr = DB2.Binary(istr)
      except:
        pass
    else:
      try:
        istr = self.unique.lookup(istr)
      except AttributeError:
        if not self.unique:
          istr = istr_orig
        else:
          raise LIGOLwParseError('unique id table has not been initialized')
    return istr

  def parsetuple(self,xmltuple):
    """
    Parse an XML tuple returned by pyRXP into a dictionary
    of LIGO metadata elements. The dictionary contains one
    entry for each table found in the XML tuple.
    """
    # first extract all the table and columns from the tuple from the
    # children of the ligo lightweight parent tuple
    table = {}
    tupleidx = 0
    for tag in xmltuple[2]:
      if tag[0] == 'Table' or tag[0] == 'table':
        tab = tag[1]['Name'].encode('ascii').lower()
        try:
          tab = self.tabrx.match(tab).group(2)
        except AttributeError:
          raise LIGOLwParseError('unable to parse a valid table name '+tab)
        # initalize the table dictionary for this table
        table[tab] = { 
          'pos' : tupleidx,
          'column' : {},
          'stream' : (), 
          'query' : ''
          }
        # parse for columns in the tables children
        # look for the column name and type in the attributes
        # store the index in which the columns were found as
        # we need this to decode the stream later
        for subtag in tag[2]:
          if subtag[0] == 'Column' or subtag[0] == 'column':
            col = subtag[1]['Name'].encode('ascii').lower()
            try:
              col = self.colrx.match(col).group(3)
            except AttributeError:
              raise LIGOLwParseError('unable to parse a valid column name '+col)
            try:
              typ = subtag[1]['Type'].encode('ascii').lower()
            except KeyError:
              raise LIGOLwParseError('type is missing for column '+col)
            table[tab]['column'][col] = typ
            table[tab].setdefault('orderedcol',[]).append(col)
      tupleidx += 1

    # now iterate the dictionary of tables we have created looking for streams
    for tab in table.keys():
      for tag in xmltuple[2][table[tab]['pos']][2]:
        if tag[0] == 'Stream' or tag[0] == 'stream':
          # store the stream delimiter and create the esacpe regex
          try:
            delim = tag[1]['Delimiter'].encode('ascii')
          except KeyError:
            raise LIGOLwParseError('stream is missing delimiter')
          if delim != ',':
            raise LIGOLwParseError('unable to handle stream delimiter: '+delim)

          # If the result set is empty tag[2] is an empty array, which causes
          # the next step to fail.  Add an empty string in this case.
          if len(tag[2]) == 0:
            tag[2].append("")


          # strip newlines from the stream and parse it
          stream = csv.reader([re.sub(r'\n','',tag[2][0])],LIGOLWStream).next()

          # turn the csv stream into a list of lists
          slen = len(stream)
          ntyp = len(table[tab]['column'])
          mlen, lft = divmod(slen,ntyp)
          if lft != 0:
            raise LIGOLwParseError('invalid stream length for given columns')
          lst = [[None] * ntyp for i in range(mlen)]

          # translate the stream data to the correct data types
          for i in range(slen):
            j, k = divmod(i,ntyp)
            try:
              thiscol = table[tab]['orderedcol'][k]
              if len( stream[i] ) == 0:
                lst[j][k] = None
              else:
                lst[j][k] = self.types[table[tab]['column'][thiscol]](stream[i])
            except (KeyError, ValueError), errmsg:
              msg = "stream translation error (%s) " % str(errmsg)
              msg += "for column %s in table %s: %s -> %s" \
                % (tab,thiscol,stream[i],str(table[tab])) 
              raise LIGOLwParseError(msg)
          table[tab]['stream'] = map(tuple,lst)

    # return the created table to the caller
    return table
          

class LIGOMetadata:
  """
  LIGO Metadata object class. Contains methods for parsing a LIGO
  lightweight XML file and inserting it into a database, executing
  and SQL query to retrive data from the database and writing it
  to a LIGO lightweight XML file
  """
  def __init__(self,xmlparser=None,lwtparser=None,ldb=None):
    """
    Connects to the database and creates a cursor. Initializes the unique
    id table for this LIGO lw document.

    ldb = LIGOMetadataDatabase object
    xmlparser = pyRXP XML to tuple parser object
    lwtparser = LIGOLwParser object (tuple parser)
    """
    self.ldb = ldb
    if self.ldb:
      self.dbcon = DB2.connect(dsn=self.ldb.database, uid='', pwd='')
      self.curs = self.dbcon.cursor()
    else:
      self.dbcon = None
      self.curs = None
    self.xmlparser = xmlparser
    self.lwtparser = lwtparser
    if lwtparser:
      self.lwtparser.unique = None
    self.table = {}
    self.strtoxml = Xlator({ r'&' : r'&amp;', r'>' : r'&gt;', r'<' : r'&lt;', '\\' : '\\\\', '\"' : '\\\"' }) # Note: see https://www.gravity.phy.syr.edu/dokuwiki/doku.php?id=rpfisher:gluebughunt if this is confusing, the parser just cleanly handles the conversion of everything

  def __del__(self):
    if self.curs:
      self.curs.close()
    if self.dbcon:
      self.dbcon.close()

  def reset(self):
    """Clear any existing table"""
    if self.table:
      del self.table
    self.table = {}

  def parse(self,xml):
    """
    Parses an XML document into a form read for insertion into the database

    xml = the xml document to be parsed
    """
    if not self.xmlparser:
      raise LIGOLwParseError("pyRXP parser not initialized")
    if not self.lwtparser:
      raise LIGOLwParseError("LIGO_LW tuple parser not initialized")
    xml = "".join([x.strip() for x in xml.split('\n')])
    ligolwtup = self.xmlparser(xml)
    if self.curs:
      self.lwtparser.unique = UniqueIds(self.curs)
    self.table = self.lwtparser.parsetuple(ligolwtup)

  def add_lfn(self,lfn):
    """
    Add an LFN table to a parsed LIGO_LW XML document.

    lfn = lfn to be added
    """
    if len(self.table['process']['stream']) > 1:
      msg = "cannot add lfn to table with more than one process"
      raise LIGOLwParseError(msg)
    # get the process_id from the process table
    pid_col = self.table['process']['orderedcol'].index('process_id')
    pid = self.table['process']['stream'][0][pid_col]
    try:
      self.table['lfn']['stream'].append((pid,lfn))
    except KeyError:
      self.table['lfn'] = {
        'pos' : 0,
        'column' : {'process_id' : 'ilwd:char', 'name' : 'lstring'},
        'stream' : [(pid, lfn)],
        'query' : '',
        'orderedcol' : ['process_id', 'name' ]
        }

  def set_dn(self,dn):
    """
    Use the domain column in the process table to store the DN

    dn = dn to be added
    """
    try:
      domain_col = self.table['process']['orderedcol'].index('domain')
      for row_idx in range(len(self.table['process']['stream'])):
        row_list = list(self.table['process']['stream'][row_idx])
        row_list[domain_col] = dn
        self.table['process']['stream'][row_idx] = tuple(row_list)
    except ValueError:
      self.table['process']['column']['domain'] = 'lstring'
      self.table['process']['orderedcol'].append('domain')
      for row_idx in range(len(self.table['process']['stream'])):
        row_list = list(self.table['process']['stream'][row_idx])
        row_list.append(dn)
        self.table['process']['stream'][row_idx] = tuple(row_list)
    
  def insert(self):
    """Insert the object into the database"""
    if not self.curs:
      raise LIGOLwDBError("Database connection not initalized")
    if len(self.table) == 0:
      raise LIGOLwDBError('attempt to insert empty table')
    for tab in self.table.keys():
      # find and add any missing unique ids
      generate = []
      missingcols = [k for k in self.ldb.uniqueids[tab] 
        if k not in self.table[tab]['column']]
      for m in missingcols:
        generate.append(',BLOB(GENERATE_UNIQUE())')
        self.table[tab]['orderedcol'].append(m)
      # and construct the sql query
      self.table[tab]['query'] = ' '.join( 
        ['INSERT INTO', tab, '(', ','.join(self.table[tab]['orderedcol']), 
        ') VALUES (', ','.join(['?' for x in self.table[tab]['column']]) , 
        ''.join(generate), ')'])
    for tabtup in self.ldb.tables:
      tab = tabtup[0].lower()
      try:
        try: 
          self.curs.executemany(self.table[tab]['query'],
            self.table[tab]['stream'])
          rowcount = self.curs.rowcount
        except DB2.Error as e:
          self.curs.execute('rollback')
          msg = e[2] 
          msg += self.xml() + '\n' 
          msg += str(self.table[tab]['query']) + '\n' 
          msg += str(self.table[tab]['stream']) + '\n'
          raise LIGOLwDBError(msg)
        except DB2.Warning, e:
          self.curs.execute('rollback')
          raise LIGOLwDBError(e[2])
        #except Exception, e:
        #  self.curs.execute('rollback')
        #  raise LIGOLwDBError, e[2]
      except KeyError:
        pass
    self.curs.execute('commit')
    return rowcount


  def select(self,sql):
    """
    Execute an SQL select statement and stuff the results into a
    dictionary.

    sql = the (case sensitve) SQL statment to execute
    """
    if not self.curs:
      raise LIGOLwDBError("Database connection not initalized")
    if len(self.table) != 0:
      raise LIGOLwDBError('attempt to fill non-empty table from database')
    ligolw = ''
    self.table = {}
    sqltypes = {
      -2 : 'ilwd:char_u',
      1 : 'lstring',
      3 : 'real_8',
      4  : 'int_4s',
      5 : 'int_2s',
      7 : 'real_4',
      8 : 'real_8',
      12 : 'lstring',
      93 : 'lstring', 
      }
    try:
      tab = re.compile(r'[Ff][Rr][Oo][Mm]\s+([A-Za-z0-0_]+)([,\s]+|$)').search(sql).group(1)
    except AttributeError:
      raise LIGOLwDBError('could not find table name in query ' + str(sql))
    self.table[tab] = {
      'pos' : 0,
      'column' : {},
      'stream' : (), 
      'query' : sql
      }
    try:
      self.curs.execute(sql)
    except DB2.Error as e:
      raise LIGOLwDBError(e[2])
    desc = self.curs.description
    for col,typ,disp,intsz,prec,sca,nul in desc:
      try:
        self.table[tab]['column'][col] = sqltypes[typ]
      except KeyError:
        raise LIGOLwDBError('unknown type returned by database ' + str(typ))
      self.table[tab].setdefault('orderedcol',[]).append(col)

    try:
      self.table[tab]['stream'] = self.curs.fetchall()
    except DB2.Error as e:
      raise LIGOLwDBError(e[2])

    return len(self.table[tab]['stream'])

  def xml(self, ilwdchar_to_hex = True):
    """Convert a table dictionary to LIGO lightweight XML"""
    if len(self.table) == 0:
      raise LIGOLwDBError('attempt to convert empty table to xml')
    ligolw = """\
<?xml version='1.0' encoding='utf-8' ?>
<?xml-stylesheet type="text/xsl" href="ligolw.xsl"?>
<!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt">
<LIGO_LW>
"""

    for tab in self.table.keys():
      try:
        ligolw += '   <Comment>'+self.strtoxml.xlat(self.table[tab]['query'])+'</Comment>\n'
      except KeyError:
        pass
      ligolw += '   <Table Name="'+tab+':table">\n'
      for col in self.table[tab]['orderedcol']:
        ligolw +='      <Column Name="'+tab.lower()+':'+col.lower()+'" Type="'+self.table[tab]['column'][col].lower()+'"/>\n'
      ligolw += '      <Stream Name="'+tab.lower()+':table" Type="Local" Delimiter=",">\n'
      stridx = 0
      ligolw += '      '
      for tup in self.table[tab]['stream']:
        if stridx != 0:
          ligolw += ',\n      '
        colidx = 0
        for tupi in tup:
          if tupi is not None:
            coltype = self.table[tab]['column'][self.table[tab]['orderedcol'][colidx]]
            if re.match(r'\Ailwd:char_u\Z',coltype):
              ligolw += '"'
              for ch in str(tupi):
                # NOTE: escape the backslash in the ilwd:char_u octal string
                ligolw += '\\\\%.3o' % (ord(ch))
              ligolw += '"'
            elif re.match(r'\Ailwd:char\Z',coltype):
              if ilwdchar_to_hex is True:
                # encode in DB2-style hex (e.g., "x'deadbeef'")
                ligolw += '"x\''
                for ch in str(tupi):
                  ligolw += "%02x" % ord(ch)
                ligolw += '\'"'
              else:
                ligolw += '"' + str(tupi) + '"'
            elif re.match(r'\Alstring\Z',coltype):
              # this santizes the contents of tupi in several ways: 
              # strtoxml.xlat escapes any double-quote and
              # backslash chars (with a preceding blackslash); and
              # then replaces <>& chars with their html
              # code equivalents
              # NOTE: string_format_func was removed so the enclosing ""
              # chars need to be added ourselves
              ligolw += '"'+self.strtoxml.xlat(tupi)+'"' 
            elif re.match(r'\Areal_4\Z',coltype):
              ligolw += '%13.7e' % tupi
            elif re.match(r'\Areal_8\Z',coltype):
              ligolw += '%22.16e' % tupi
            else:
              ligolw += str(tupi)
          else:
            ligolw += ''
          if colidx < (len(self.table[tab]['column']) - 1):
            ligolw += ','
          colidx += 1
        stridx += 1
      ligolw += '\n      </Stream>\n'
      ligolw += '   </Table>\n'
    ligolw += '</LIGO_LW>'

    return ligolw
