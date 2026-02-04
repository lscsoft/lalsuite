# Author: Karl Wette, 2026

"""LALSuite library build configuration interface"""

import re
from pathlib import Path


class LALLibraryBuildConfig:
    """LALSuite library build configuration"""

    name: str
    version: str

    @staticmethod
    def _split_args(s):

        # Split an Autoconf function argument list
        tt = [""]
        d = 0
        for c in s:
            if c == ",":
                if d == 0:
                    tt.append("")
                else:
                    tt[-1] += c
            elif c == "[":
                d += 1
            elif c == "]":
                d -= 1
            else:
                tt[-1] += c
        return [t.strip() for t in tt]

    def __init__(self, lallib, rootdir=Path.cwd()):

        # Read Autoconf configure.ac
        if lallib == "lalsuite":
            self._cfg_ac_file = rootdir / "configure.ac"
        else:
            self._cfg_ac_file = rootdir / lallib / "configure.ac"
        self._ac_lines = []
        self._ac_subst = {}
        ac_vars = {}
        with self._cfg_ac_file.open("r") as f:
            for line in f:
                line = line.rstrip()

                # Read AC_INIT line
                m = re.fullmatch(r"AC_INIT\((.*)\)", line)
                if m is not None:
                    args = self._split_args(m.group(1))
                    self.name, self.version, self._email = args
                    self._ac_lines.append(("AC_INIT", None))
                    continue

                # Read AC_SUBST line
                m = re.fullmatch(r"AC_SUBST\((.*)\)", line)
                if m is not None:
                    args = self._split_args(m.group(1))
                    var = args.pop(0)
                    try:
                        val = args.pop(0)
                    except IndexError:
                        val = None
                    self._ac_subst[var] = val
                    self._ac_lines.append(("AC_SUBST", var))
                    continue

                # Read variable assignments
                m = re.fullmatch(r"([A-Z_]+)=(.*)", line)
                if m is not None:
                    var = m.group(1)
                    val = m.group(2)
                    ac_vars[var] = val
                    self._ac_lines.append(("var", (var, val)))
                    continue

                # Read all other lines
                self._ac_lines.append(("line", line))

        # Set AC_SUBST variable values
        for var, val in ac_vars.items():
            if var in self._ac_subst:
                self._ac_subst[var] = val

    def __contains__(self, var):
        return var in self._ac_subst

    def __getitem__(self, var):
        return self._ac_subst[var]

    def __setitem__(self, var, val):
        self._ac_subst[var] = str(val)

    def write(self):

        # Write Autoconf configure.ac
        with self._cfg_ac_file.open("w") as f:
            for key, args in self._ac_lines:

                # Write AC_INIT line
                if key == "AC_INIT":
                    print(
                        f"AC_INIT([{self.name}],[{self.version}],[{self._email}])",
                        file=f,
                    )
                    continue

                # Write AC_SUBST line
                if key == "AC_SUBST":
                    var = args
                    val = self._ac_subst[var]
                    if val is None:
                        print(f"AC_SUBST([{var}])", file=f)
                    else:
                        print(f"AC_SUBST([{var}],[{val}])", file=f)
                    continue

                # Write variable assignments, if not already written by AC_SUBST
                if key == "var":
                    (var, val) = args
                    if var not in self._ac_subst:
                        print(f"{var}={val}", file=f)
                    continue

                # Write all other lines
                line = args
                print(line, file=f)
