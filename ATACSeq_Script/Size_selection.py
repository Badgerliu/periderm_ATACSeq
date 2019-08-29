#!/usr/bin/env python
try:
    import pysam
    import sys
except ImportError:
    print("make sure you have pysam and sys modules installed")

bamfile = sys.argv[1]
samfile = pysam.Samfile(bamfile, 'rb')

## limits, from Buenrosto 2013, 'Nucleosome positioning.'
nucfree = range(0, 100)
nucfreefile = 'nucFreeReads.' + bamfile + '.sam'
nucfreesam = pysam.Samfile(nucfreefile, "wh", header=samfile.header)

mononuc = range(180, 248)
mononucfile = 'monoNucReads.' + bamfile + '.sam'
mononucsam = pysam.Samfile(mononucfile, "wh", header=samfile.header)


print('Splitting reads by TLEN fields in ' + bamfile + '...')

for alignedread in samfile.fetch():
  if abs(alignedread.tlen) in nucfree: 
    nucfreesam.write(alignedread)
  elif abs(alignedread.tlen) in mononuc: 
    mononucsam.write(alignedread)

nucfreesam.close()
mononucsam.close()
