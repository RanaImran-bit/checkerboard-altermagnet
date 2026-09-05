# Give every run its own random seed.
#
# cp.f90 sets  ISEED = -valu(8)*(cpuRank+1),  and valu(8) is the MILLISECOND
# field of date_and_time. At NP=8 the cpuRank factor spreads the ranks out, but
# the campaign runs one rank per core, so cpuRank is always 0 and the seed can
# take at most 1000 values. Launching hundreds of runs at once then hands the
# same seed to different runs, and averaging over "independent" seeds that are
# actually identical understates the error bars. That is a silent statistical
# error, not a crash.
#
# With this patch the launcher writes a unique integer into seed.dat and the code
# uses it. The clock path stays as the fallback so the file is optional.
import sys
s = open("cp.f90").read()
old = """      ISEED=(-valu(8)*(cpuRank+1))
"""
new = """      ISEED=(-valu(8)*(cpuRank+1))
!! CAMPAIGN: prefer a unique per-run seed from seed.dat. The clock-derived value
!! above spans only the millisecond field, so concurrent one-rank runs collide.
      open(77,file='seed.dat',status='old',iostat=iosSeed)
      if (iosSeed == 0) then
         read(77,*,iostat=iosSeed) iseedFile
         close(77)
         if (iosSeed == 0) ISEED = -abs(iseedFile) - cpuRank
      end if
      if (ISEED >= 0) ISEED = -12345 - cpuRank
      write(*,*) 'ISEED=',ISEED
      flush(6)
"""
if s.count(old) != 1:
    sys.exit(f"expected 1 occurrence of the ISEED line, found {s.count(old)}")
open("cp.f90","w").write(s.replace(old, new))

p = open("parameter.f90").read()
anchor = "integer(k4b)::nseed(4),ISEED"
if p.count(anchor) != 1:
    sys.exit("parameter.f90 anchor not found")
open("parameter.f90","w").write(p.replace(anchor,
    anchor + "\ninteger::iosSeed,iseedFile      ! per-run seed read from seed.dat"))
print("seed fix applied to cp.f90 and parameter.f90")
