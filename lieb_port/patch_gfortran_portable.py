# Make mc2duph.f90 build under gfortran as well as ifort.
#
# The three kNames tables are written as (/ 'gx_up', 'gx_dn', ... /) with literals
# of differing length. The standard requires every element of a character array
# constructor to have the SAME length; ifort silently pads, gfortran rejects it:
#
#   Error: Different CHARACTER lengths (5/4) in array constructor
#
# and then every trim(kNames(j)) downstream fails too, because the failed
# declaration leaves kNames implicitly typed real.
#
# The fix is the F2003 type-spec form, [character(len=32) :: ...], which pads
# automatically and is accepted by both compilers. Nothing about the values or
# the file names changes.
#
# Needed because 113new has no Intel oneAPI for this account and builds with
# GNU Fortran + OpenMPI.
import re, sys

src = open("mc2duph.f90").read()
pat = re.compile(r"(character\(len=(\d+)\), parameter :: kNames\(\d+\) = )\(/(.*?)/\)", re.S)
n = 0
def fix(m):
    global n; n += 1
    return f"{m.group(1)}[character(len={m.group(2)}) ::{m.group(3)}]"
out = pat.sub(fix, src)
if n != 3:
    sys.exit(f"expected 3 kNames constructors, rewrote {n}")
open("mc2duph.f90", "w").write(out)
print(f"rewrote {n} kNames constructors to the type-spec form")
