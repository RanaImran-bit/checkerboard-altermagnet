# Make the runs observable while they run.
#
# Neither progress channel is flushed, so both look dead until the job exits:
#   cp.f90       writes 'equilibrate: N' to unit 5, which is time.dat
#   mc2duph.f90  writes the per-block Egrowth / W_Up / W_Down line to stdout
# A block line is about 60 bytes and there are only 10 of each, so an 8 KB
# Fortran buffer never fills and nothing reaches disk. That is why a job 46 hours
# into its walk showed a 0-byte time.dat and no block lines, and why I misread a
# working run as a stalled one.
import sys

def sub1(path, old, new, tag):
    s = open(path).read()
    if s.count(old) != 1:
        sys.exit(f"[{tag}] expected 1 occurrence in {path}, found {s.count(old)}")
    open(path, "w").write(s.replace(old, new))

sub1("cp.f90",
     "        write(5,*) 'equilibrate:',iblk\n",
     "        write(5,*) 'equilibrate:',iblk\n"
     "        flush(5)\n"
     "        if (myID==0) then\n"
     "          write(*,*) 'equilibrate block',iblk,'of',nblkeq\n"
     "          flush(6)\n"
     "        end if\n",
     "equilibrate")

sub1("mc2duph.f90",
     "write(*,*) iblk,eavg,sum(wgtwlkr)\n",
     "write(*,*) iblk,eavg,sum(wgtwlkr)\nflush(6)\n",
     "block line")
print("flush added: equilibration progress and per-block line now appear live")
