# generate_swig_iface.py - generate SWIG interface
#
# Copyright (C) 2011--2014 Karl Wette
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with with program; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

__author__ = 'Karl Wette <karl.wette@ligo.org>'
__copyright__ = 'Copyright (C) 2011--2014 Karl Wette'

import sys, os, re
from xml.etree.cElementTree import ElementTree

# print error message and exit
def fail(msg):
    sys.stderr.write('%s: %s\n' % (sys.argv[0], msg))
    sys.exit(1)

# get input arguments
package_name, dependencies, symbol_prefixes, input_xml, output_iface = sys.argv[1:]
assert(package_name[:3] == 'LAL')

# parse XML input
tree = ElementTree()
try:
    tree.parse(input_xml)
except:
    fail("could not parse XML input from '%s'" % input_xml)

# function: return value corresponding to 'name' in SWIG XML <attributelist>
def get_swig_attr(elem, name):
    if elem == None:
        return None
    for attr_elem in elem.findall('attributelist/attribute'):
        if attr_elem.get("name", "") == name:
            return attr_elem.get('value')
    return None

# find root of preprocessing interface parse tree
for incl_elem in tree.findall('include'):
    if get_swig_attr(incl_elem.find('module'), 'name') == 'swiglal_preproc':
        root = incl_elem
        break
else:
    fail('could not find root of preprocessing interface parse tree')

# get list of interface headers
headers = dict()
for incl_elem in root.findall('include'):

    # find header filename
    header_name = get_swig_attr(incl_elem, 'name')
    if header_name == None:
        fail('could not find filename of header file')
    _, header_name = os.path.split(header_name)

    # filter out SWIG interface headers
    if header_name.endswith('.i'):
        continue

    # add this header
    headers[header_name] = incl_elem

# get dictionaries of symbols from interface headers
functions = dict()
funcnames = dict()
tdstructs = dict()
tagnames = dict()
structs = dict()
clear_macros = dict()
for header_name in headers:
    clear_macros[header_name] = set()

    # C declarations
    for cdecl in headers[header_name].findall('cdecl'):
        cdecl_name = get_swig_attr(cdecl, 'name')
        cdecl_kind = get_swig_attr(cdecl, 'kind')
        cdecl_type = get_swig_attr(cdecl, 'type')
        if cdecl_name == None:
            fail("cdecl in header '%s' has no name" % header_name)
        if cdecl_kind == None:
            fail("cdecl '%s' in header '%s' has no kind" % (cdecl_name, header_name))
        if cdecl_type == None or cdecl_type == '':
            fail("cdecl '%s' in header '%s' has no type" % (cdecl_name, header_name))

        # SWIGLAL() and SWIGLAL_CLEAR() macros
        if cdecl_name in ['__swiglal__', '__swiglal_clear__']:
            cdecl_value = get_swig_attr(cdecl, 'value')
            if cdecl_value == None:
                fail("cdecl '%s' in header '%s' has no value" % (cdecl_name, header_name))
            macro = re.sub(r'\s', '', cdecl_value)
            if cdecl_name == '__swiglal__':
                if macro in clear_macros[header_name]:
                    fail("duplicate definition of SWIGLAL(%s)", macro)
                else:
                    clear_macros[header_name].add(macro)
            else:
                if macro in clear_macros[header_name]:
                    clear_macros[header_name].remove(macro)
                else:
                    fail("cannot clear undefined macro SWIGLAL(%s)", macro)

        # functions
        elif cdecl_kind == 'function':
            if cdecl_name in functions:
                fail("duplicate function '%s' in header '%s'" % (cdecl_name, header_name))
            functions[cdecl_name] = cdecl
            funcnames[cdecl_name] = cdecl_name

        # typedefs to structs
        elif cdecl_kind == 'typedef' and cdecl_type.startswith('struct '):
            if cdecl_name in tdstructs:
                fail("duplicate struct typedef '%s' in header '%s'" % (cdecl_name, header_name))
            tagnames[cdecl_name] = cdecl_type[7:]
            tdstructs[cdecl_name] = cdecl

    # structs
    for struct in headers[header_name].findall('class'):
        struct_name = get_swig_attr(struct, 'name')
        if struct_name == None:
            fail("struct in header '%s' has no name" % header_name)
        if get_swig_attr(struct, 'unnamed') != None:
            fail("struct '%s' in header '%s' has no tag-name" % (struct_name, header_name))
        if struct_name in structs:
            fail("duplicate struct '%s' in header %s" % (struct_name, header_name))
        structs[struct_name] = struct

# list of symbol prefixes, for renaming
symbol_prefix_list = symbol_prefixes.split()
symbol_prefix_list.append('')

# function: build renaming map for symbols, using symbol prefix order
def rename_symbols(rename_taken, rename_map, symbols, name_map):

    # rank symbol by where their prefix comes in prefix list
    rename_rank = dict()
    for symbol_key in symbols:
        for rank in range(len(symbol_prefix_list)):
            if symbol_key.startswith(symbol_prefix_list[rank]):
                rename_rank[symbol_key] = rank
                break

    # iterate over symbols in order of decreasing rank
    for symbol_key in sorted(symbols, key = lambda k : rename_rank[k]):

        # get name of symbol
        symbol_name = name_map[symbol_key]

        # check for duplicate symbols
        if symbol_name in rename_map:
            fail("duplicate symbols '%s' in interface" % symbol_key)

        # ignore symbol if the 'feature_ignore' attribute is set,
        # otherwise strip prefix from symbol key to get re-name
        if get_swig_attr(symbols[symbol_key], 'feature_ignore') != None:
            symbol_rename = None
        else:
            rank = rename_rank[symbol_key]
            symbol_rename = symbol_key[len(symbol_prefix_list[rank]):]

        # ignore symbol if another symbol has taken the re-name,
        # otherwise record that a symbol has taken this re-name
        if symbol_rename != None:
            if symbol_rename in rename_taken:
                symbol_rename = None
            else:
                rename_taken.add(symbol_rename)

        # if re-name differs from original name, add to rename map
        if symbol_rename != symbol_name:
            rename_map[symbol_name] = symbol_rename

# rename functions and struct typedefs, using symbol prefix order
rename_taken = set()
rename_map = dict()
rename_symbols(rename_taken, rename_map, functions, funcnames)
rename_symbols(rename_taken, rename_map, tdstructs, tagnames)

# look for a destructor function for each struct
dtor_functions = dict()
dtor_name_regexp = re.compile(r'(Destroy|Free|Close)([A-Z0-9_]|$)')
dtor_decl_regexp = re.compile(r'^f\(p\.(.*)\)\.$')
for function_name in functions:

    # function must match destructor name regexp
    dtor_name_match = dtor_name_regexp.search(function_name)
    if dtor_name_match == None:
        continue

    # function must take a single pointer argument
    dtor_decl = get_swig_attr(functions[function_name], 'decl')
    dtor_decl_match = dtor_decl_regexp.match(dtor_decl)
    if dtor_decl_match == None:
        continue
    dtor_struct_name = dtor_decl_match.group(1)

    # function argument must be a struct name
    if not dtor_struct_name in tdstructs:
        continue

    # function must return either void or int; fail otherwise
    dtor_retn_type = get_swig_attr(functions[function_name], 'type')
    if not dtor_retn_type in ['void', 'int']:
        fail("destructor function '%s' has invalid return type '%s'" % (function_name, dtor_retn_type))

    # struct must not already have a destructor; fail otherwise
    if dtor_struct_name in dtor_functions:
        fail("struct typedef '%s' has duplicate destructors '%s' and '%s'" % (dtor_struct_name, function_name, dtor_functions[dtor_struct_name]))

    # save destructor name
    dtor_functions[dtor_struct_name] = function_name

    # remove destructor function from interface
    rename_map[function_name] = None

# start SWIG interface file
try:
    f = open(output_iface, 'w')
except:
    fail("could not open output interface file '%s'" % output_iface)
f.write('// %s: generated by %s\n' % (output_iface, sys.argv[0]))
f.write('%%module %s;\n' % package_name.lower())
f.write('%include <lal/SWIGCommon.i>\n')

# import dependent modules
f.write('#ifndef SWIGIMPORTED\n')
for module in dependencies.split():
    f.write('%%import <lal/swig_%s.i>\n' % module.lower())
f.write('#endif // !SWIGIMPORTED\n')

# include interface headers in wrapping code
f.write('%header %{\n')
for header_name in sorted(headers):
    f.write('#include <lal/%s>\n' % header_name)
f.write('%}\n')

# perform symbol renames
for name in sorted(rename_map):
    rename = rename_map[name]
    if rename == None:
        f.write('%%ignore %s;\n' % name)
    else:
        f.write('%%rename("%s") %s;\n' % (rename, name))

# perform operations on functions
func_retn_1starg_regexp = re.compile(r'^f\(((?:p\.)+[^,]+?)(?:,[^,]+)*\)\.\1$')
for function_name in sorted(functions):

    # skip ignored functions
    if rename_map.get(function_name, '') == None:
        continue

    # apply %newobject to all functions
    f.write('%%newobject %s;\n' % function_name)

    # indicate if the return type of a function is a pointer type, and matches the
    # type of its first argument; many LAL functions return their first argument
    # after performing some operation on it, but we want to prevent two different
    # SWIG wrapping objects from owning the same LAL memory
    func_decl = get_swig_attr(functions[function_name], 'decl')
    func_type = get_swig_attr(functions[function_name], 'type')
    func_decl_type = func_decl + func_type
    if func_retn_1starg_regexp.match(func_decl_type) != None:
        f.write('%%header %%{#define swiglal_return_1starg_%s%%}\n' % function_name)

# include interface headers, and clear SWIGLAL() macros afterwards
f.write('%%include <lal/SWIG%sAlpha.i>\n' % package_name)
for header_name in sorted(headers):
    f.write('%%include <lal/%s>\n' % header_name)
    for macro_name in clear_macros[header_name]:
        f.write('SWIGLAL_CLEAR(%s);\n' % macro_name)
f.write('%%include <lal/SWIG%sOmega.i>\n' % package_name)

# create constructors and destructors for structs
for struct_name in sorted(tdstructs):
    struct_tagname = tagnames[struct_name]
    struct_opaque = '%i' % (not struct_tagname in structs)
    struct_dtor_function = dtor_functions.get(struct_name, '')
    struct_args = (struct_name, struct_tagname, struct_opaque, struct_dtor_function)
    f.write('%%swiglal_struct_create_cdtors(%s)\n' % ','.join(struct_args))
