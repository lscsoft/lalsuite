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
# Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301  USA

"""Generate the swig<package>.i interface file for a given LAL subpackage

This script scans the SWIG pre-processed XML definitions for the given library
and maps each symbol to a SWIG definition.
This also handles on-the-fly symbol renaming, mainly to make the SWIG interface
more natural for high-level languages that support namespaces.
"""

from __future__ import print_function

__author__ = "Karl Wette <karl.wette@ligo.org>"
__copyright__ = 'Copyright (C) 2011--2014 Karl Wette'

import argparse
import os
import re
import sys
from xml.etree.cElementTree import ElementTree
if sys.version_info >= (3, 7):
    OrderedDict = dict
else:
    from collections import OrderedDict

# renaming prefixes: these are stripped from the object names
#                    provided in the bindings
CONSTANT_RENAMES = (
    "XLAL_",
    "LAL_",
)
VARIABLE_RENAMES = (
    "lal",
)

# functions that return their first input are modified
FUNC_RETN_1STARG_REGEX = re.compile(r'^f\(((?:p\.)+[^,]+?)(?:,[^,]+)*\)\.\1$')

# destructors functions are linked as __del__ methods of the relevant
# classes, and not exposed individually in the bindings
DESTRUCTOR_REGEX = re.compile(r'(Destroy|Free|Close)([A-Z0-9_]|$)')
DESTRUCTOR_DECL_REGEX = re.compile(r'^f\(p\.(.*)\)\.$')

# old LAL error messages and codes are ignored from the bindings entirely
OLD_LAL_EMSG_REGEX = re.compile(r'([A-Z0-9_]+H_)MSG(E[A-Z0-9]+)$')


# -- symbol parsing -----------------------------

class SwigSymbols(object):
    """Object to represent the bundle of symbols parsed from a library
    to be exposed through a swig binding
    """
    def __init__(self):
        # parsing
        self.clear_macros = {}
        self.constants = {}
        self.functions = {}
        self.structs = {}
        self.tdstructs = {}
        self.tdstruct_names = {}
        self.typedefs = {}
        self.variables = {}

        # renames
        self.renames = {}
        self.rename_kind = {}

        # destructors
        self.destructors = {}

    def parse_headers(self, headers):
        """Parse all headers
        """
        for header in headers:
            self.parse_header(header, headers[header])

    def parse_header(self, header, elem):
        """Parser a single header from an XML element
        """
        self.clear_macros[header] = set()

        # declarations
        for cdecl in elem.findall('cdecl'):
            self._parse_cdecl(cdecl, header)

        # structs
        for struct in elem.findall('class'):
            self._parse_struct(struct)

        # preprocessor constants
        for constant in elem.findall('constant'):
            self._parse_constant(constant)

        # enumerations and enumeration constants
        for enum in elem.findall('enum'):
            self._parse_enum(enum)

    # -- element parsers ----

    def _parse_cdecl(self, cdecl, header):
        """Parse a C declaration
        """
        cdecl_name = get_swig_attr(cdecl, 'name')
        cdecl_kind = get_swig_attr(cdecl, 'kind')
        cdecl_type = get_swig_attr(cdecl, 'type')
        if not cdecl_name:
            raise ValueError(
                "cdecl in header '{}' has no name".format(header),
            )
        if not cdecl_kind:
            raise ValueError(
                "cdecl '{}' in header '{}' has no kind".format(
                    cdecl_name,
                    header,
                ),
            )
        if not cdecl_type:
            raise ValueError(
                "cdecl '{}' in header '{}' has no type".format(
                    cdecl_name,
                    header,
                ),
            )

        # SWIGLAL() and SWIGLAL_CLEAR() macros
        if cdecl_name in ['__swiglal__', '__swiglal_clear__']:
            self._parse_swiglal_macro(cdecl, self.clear_macros[header])

        # functions
        elif cdecl_kind == 'function':
            self._parse_function(cdecl)

        # typdefs
        elif cdecl_kind == "typedef":
            self._parse_typedef(cdecl)

        # variables
        elif (
                cdecl_kind == "variable" and
                cdecl_type not in {"SWIGLAL", "SWIGLAL_CLEAR"}
        ):
            self._parse_variable(cdecl)

    def _parse_constant(self, constant):
        """Parse a constant
        """
        self.constants[get_swig_attr(constant, 'name')] = constant

    def _parse_enum(self, enum):
        """Parse an enum (special case of constant)
        """
        if get_swig_attr(enum, 'unnamed') is not None:
            ename = get_swig_attr(enum, 'name')
            raise ValueError("enum '{}' has no tag-name".format(ename))
        for item in enum.findall('enumitem'):
            self._parse_constant(item)

    def _parse_function(self, cdecl):
        """Parse a function
        """
        cdecl_name = get_swig_attr(cdecl, 'name')
        if cdecl_name in self.functions:
            raise ValueError("duplicate function '{}'".format(cdecl_name))
        self.functions[cdecl_name] = cdecl

    def _parse_struct(self, struct):
        """Parse a struct
        """
        name = get_swig_attr(struct, 'name')
        if not name:
            raise ValueError("struct has no name")
        if get_swig_attr(struct, 'unnamed') is not None:
            raise ValueError(
                "struct '{}' has no tag-name".format(name)
            )
        if name in self.structs:
            raise ValueError("duplicate struct '{}'".format(name))
        self.structs[name] = struct

    def _parse_swiglal_macro(self, cdecl, clear_macros):
        """Parse a swiglal macro
        """
        cdecl_name = get_swig_attr(cdecl, 'name')
        cdecl_value = get_swig_attr(cdecl, 'value')
        if not cdecl_value:
            raise ValueError("cdecl '{}' has no value".format(cdecl_name))
        macro = re.sub(r'\s', '', cdecl_value)
        if cdecl_name == '__swiglal__':
            if macro in clear_macros:
                raise ValueError(
                    "duplicate definition of SWIGLAL({})".format(macro),
                )
            clear_macros.add(macro)
        else:
            if macro not in clear_macros:
                raise ValueError(
                    "cannot clear undefined macro SWIGLAL({})".format(macro),
                )
            clear_macros.remove(macro)

    def _parse_typedef(self, cdecl):
        """Parse a ``typedef``
        """
        cdecl_name = get_swig_attr(cdecl, 'name')
        cdecl_type = get_swig_attr(cdecl, 'type')
        if cdecl_type.startswith('struct '):
            if cdecl_name in self.tdstructs:
                raise ValueError(
                    "duplicate struct typedef '{}'".format(cdecl_name),
                )
            self.tdstruct_names[cdecl_name] = cdecl_type[7:]
            self.tdstructs[cdecl_name] = cdecl

        # other typedefs
        else:
            self.typedefs[cdecl_name] = cdecl_type

    def _parse_variable(self, cdecl):
        """Parse a variable
        """
        self.variables[get_swig_attr(cdecl, 'name')] = cdecl

    # -- renaming -----------

    def rename(self, prefixes):
        """Scan the parsed symbols and generating a renaming map
        according to some pre-defined rules
        """
        self._rename_old_lal_error_codes_and_messages(
            self.constants,
        )
        self._rename_symbols(
            "constant",
            self.constants,
            CONSTANT_RENAMES,
        )
        self._rename_symbols(
            "variable",
            self.variables,
            VARIABLE_RENAMES,
        )
        self._rename_symbols(
            "function",
            self.functions,
            prefixes,
        )
        self._rename_symbols(
            "tdstruct",
            self.tdstructs,
            prefixes,
            names=self.tdstruct_names,
        )
        return self.renames, self.rename_kind

    def _rename_old_lal_error_codes_and_messages(self, constants):
        """Mark old-style LAL error codes and message to be ignored
        """
        for oldemsg in self.constants:
            match = OLD_LAL_EMSG_REGEX.match(oldemsg)
            if not match:
                continue
            oldeval = ''.join(match.groups())
            if oldeval in self.constants:
                self.renames[oldemsg] = None
                self.renames[oldeval] = None

    def _rename_symbols(self, kind, symbols, prefixes, names=None):
        """Rename a class of symbols
        """
        nprefixes = len(prefixes)

        # rank symbol by where their prefix comes in prefix list
        rank = dict()
        for key in symbols:
            for i, prefix in enumerate(prefixes):
                if key.startswith(prefix):
                    rank[key] = i
                    break
            else:
                rank[key] = nprefixes

        # iterate over symbols in order of decreasing rank
        for key in sorted(symbols, key=lambda k: rank[k]):
            # ignore symbol if the 'feature_ignore' attribute is set
            if get_swig_attr(symbols[key], 'feature_ignore') is not None:
                continue

            # get name of symbol
            name = names[key] if names else key

            # check for duplicate symbol renames
            if self.renames.get(name):
                raise ValueError(
                    "duplicate symbols '%s' in interface".format(key),
                )

            # strip prefix from symbol key to get re-name,
            # so long as re-name would not start with a digit or '_'
            rename = key
            if rank[key] < nprefixes:
                prefix_len = len(prefixes[rank[key]])
                if key[prefix_len] not in "0123456789_":
                    rename = key[prefix_len:]

            # if re-name has already been taken
            if rename in self.rename_kind:
                # if the re-name has been taken by a symbol of the same kind,
                # then ignore this symbol; otherwise, do not rename symbol
                if self.rename_kind[rename] == kind:
                    rename = None
                else:
                    rename = name

            # otherwise make this symbol re-name as having been taken
            else:
                self.rename_kind[rename] = kind

            # if re-name differs from original name, add to rename map
            if rename != name:
                self.renames[name] = rename

    # dectructors
    def find_destructors(self):
        """Find and map all object destructors
        """
        for name in filter(DESTRUCTOR_REGEX.search, self.functions):
            self._find_destructor(name, self.functions[name])
        return self.destructors

    def _find_destructor(self, name, elem):
        """Find and map a desctructor for a single object
        """
        # function must take a single pointer argument
        decl = get_swig_attr(elem, 'decl')
        decl_match = DESTRUCTOR_DECL_REGEX.match(decl)
        if not decl_match:
            return
        struct_name = decl_match.group(1)

        # resolve any nested typedefs
        while struct_name in self.typedefs:
            struct_name = self.typedefs[struct_name]

        # function argument must be a struct name
        if struct_name not in self.tdstructs:
            return

        # function must return either void or int; fail otherwise
        returntype = get_swig_attr(elem, 'type')
        if returntype not in {'void', 'int'}:
            raise ValueError(
                "destructor function '{}' has invalid return "
                "type '{}'".format(name, returntype),
            )

        # if struct already has a destructor, prefer the destructor
        # containing the struct name, otherwise fail
        if (
                struct_name in self.destructors and
                struct_name in self.destructors[struct_name]
        ):
            if struct_name not in name:
                return
            raise ValueError(
                "struct typedef '{}' has duplicate destructors '{}' and "
                "'{}'".format(
                    struct_name,
                    name,
                    self.destructors[struct_name],
                ),
            )

        # save destructor name
        self.destructors[struct_name] = name

        # remove destructor function from interface
        self.renames[name] = None


# -- utilities ----------------------------------

def parse_input_xml(path):
    """Parse a SWIG pre-processed XML file
    """
    tree = ElementTree()
    try:
        tree.parse(path)
    except SyntaxError as exc:
        exc.args = (
            "could not parse XML input from '{}': {}".format(str(exc), path),
        )
        raise

    # find root of preprocessing interface parse tree
    for incl in tree.findall('include'):
        if get_swig_attr(incl.find('module'), 'name') == 'swiglal_preproc':
            return incl
    raise RuntimeError(
        "could not find root of preprocessing interface parse tree",
    )


def get_header_list(elem):
    """Get list of interface headers, ordered by #include precedence
    """
    headers = OrderedDict()
    for incl in elem.findall('include'):
        # find header filename
        header_name = get_swig_attr(incl, 'name')
        if header_name is None:
            raise RuntimeError("could not find filename of header file")
        header_name = os.path.basename(header_name)

        # filter out SWIG interface headers
        if header_name.endswith('.i'):
            continue

        # skip in already in header list
        if header_name in headers:
            continue

        # add headers #included by this header to list
        headers.update(get_header_list(incl))

        # add this header to list
        headers[header_name] = incl

    return headers


def get_swig_attr(elem, name):
    """Return value corresponding to 'name' in SWIG XML <attributelist>
    """
    if elem is None:
        return None
    for attr_elem in elem.findall('attributelist/attribute'):
        if attr_elem.get("name", "") == name:
            return attr_elem.get('value')


def write_interface_file(file, package, symbols, dependencies, headers):
    """Write SWIG interface file
    """
    def fprint(*args, **kwargs):
        kwargs["file"] = file
        return print(*args, **kwargs)

    fprint('// {}: generated by {}'.format(file.name, sys.argv[0]))
    fprint('%module {};'.format(package.lower()))
    fprint('%include <lal/SWIGCommon.i>')

    # import dependent modules
    fprint('#ifndef SWIGIMPORTED')
    for module in dependencies or []:
        fprint('%import <lal/swig{}.i>'.format(module.lower()))
    fprint('#endif // !SWIGIMPORTED')

    # include interface headers in wrapping code
    fprint('%header %{')
    for header in headers:
        fprint('#include <lal/%s>' % header)
    fprint('%}')

    # add renames
    for name in sorted(symbols.renames):
        if symbols.renames[name] is None:
            fprint('%ignore {};'.format(name))
    for kind in sorted(set(symbols.rename_kind.values())):
        fprint('#ifdef SWIGLAL_MODULE_RENAME_{}S'.format(kind.upper()))
        for name in sorted(symbols.renames):
            rename = symbols.renames[name]
            if rename is not None and symbols.rename_kind[rename] == kind:
                fprint('%rename("{}") {};'.format(rename, name))
        fprint('#endif // SWIGLAL_MODULE_RENAME_{}S'.format(kind.upper()))

    # perform operations on functions
    for function_name in sorted(symbols.functions):
        # skip ignored functions
        if symbols.renames.get(function_name, '') is None:
            continue

        # apply %newobject to all functions
        fprint('%newobject {};'.format(function_name))

        # get function declaration and type
        func_decl = get_swig_attr(symbols.functions[function_name], 'decl')
        func_type = get_swig_attr(symbols.functions[function_name], 'type')
        func_decl_type = func_decl + func_type

        # apply the SWIGLAL_MAYBE_RETURN_INT typemap to all functions which
        # return 'int'; this causes the 'int' return to be ignored in most
        # cases, unless the typemap SWIGLAL(RETURN_VALUE(int, ...)) is
        # applied to disable this behaviour
        if func_type == 'int':
            fprint(
                '%apply int SWIGLAL_MAYBE_RETURN_INT {{ int {} }};'.format(
                    function_name,
                ),
            )

        # indicate if the return type of a function is a pointer type, and
        # matches the type of its first argument; many LAL functions return
        # their first argument after performing some operation on it, but we
        # want to prevent two different SWIG wrapping objects from owning the
        # same LAL memory
        if FUNC_RETN_1STARG_REGEX.match(func_decl_type):
            fprint(
                '%header %{{#define swiglal_return_1starg_{}%}}'.format(
                    function_name,
                ),
            )

    # include interface headers, and clear SWIGLAL() macros afterwards
    fprint('%include <lal/SWIG{}Alpha.i>'.format(package))
    for header in headers:
        fprint('%include <lal/{}>'.format(header))
        for macro_name in symbols.clear_macros[header]:
            fprint('SWIGLAL_CLEAR({});'.format(macro_name))
    fprint('%include <lal/SWIG{}Omega.i>'.format(package))

    # add extensions to structs, e.g. constructors and destructors
    for struct_name in sorted(symbols.tdstructs):
        struct_tagname = symbols.tdstruct_names[struct_name]
        struct_opaque = str(int(struct_tagname not in symbols.structs))
        struct_dtor_function = symbols.destructors.get(struct_name, '')
        fprint(
            '%swiglal_struct_extend({},{},{})'.format(
                struct_tagname,
                struct_opaque,
                struct_dtor_function,
            ),
        )


# -- parse command line -------------------------

def create_parser():
    """Create an `argparse.ArgumentParser` for this script
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument("package_name", help="name of package")
    parser.add_argument("dependencies", help="list of package dependencies")
    parser.add_argument("function_prefixes", help="function prefices")
    parser.add_argument("input_xml", help="preprocessed XML file from swig")
    parser.add_argument("output_iface", help="target output filename")
    return parser


def parse_arguments():
    """Parse command line arguments
    """
    parser = create_parser()
    args = parser.parse_args()
    assert args.package_name.startswith("LAL")
    args.dependencies = args.dependencies.split()
    args.function_prefixes = args.function_prefixes.split()
    return args


# -- run ----------------------------------------

if __name__ == "__main__":
    args = parse_arguments()

    # parse input XML file
    root = parse_input_xml(args.input_xml)

    # get list of interface headers
    headers = get_header_list(root)

    # parse all headers
    symbols = SwigSymbols()
    symbols.parse_headers(headers)

    # rename symbols
    symbols.rename(args.function_prefixes)

    # find destructors
    symbols.find_destructors()

    # write interface file
    with open(args.output_iface, "w") as file:
        write_interface_file(
            file,
            args.package_name,
            symbols,
            args.dependencies,
            headers,
        )
