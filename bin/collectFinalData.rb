#!/usr/bin/env ruby
# Copyright (C) 2007 Matt Hibbs
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.

# This is the Troilkatt version as of June 15, 2011 with changes by peak@princeton.edu including:
#  - modifications to allow use with UTF-8 input using either Ruby 1.8 or 1.9
#  - -Xmx* option for java; ARGV[2] as "="; / in paths is now optional

BN = File.basename($0)

@@verbosity = false

jopt="-Xmx1g"

while ARGV[0] =~ /^-/
  case ARGV[0]
  when "-v", "--verbose"
    @@verbosity = true
    ARGV.shift
  when /^-Xmx/
    jopt=ARGV[0]
    ARGV.shift
  end
end


if ARGV.length < 4
  print <<-EOH
Syntax: #{BN} [OPTIONS] ARGS
ARGS:
   0 - input pcl pathname or filename; the file name must have a period-delimited prefix
   1 - info file (aka configuration file - often GSEnnn.info or allInfo.txt)
   2 - path to DivLogNorm.jar, or = if the .jar file is in the same directory as this script
   3 - path to final data folder, referred to as OUTDIR below

Options:
   -v | --verbose :: verbose (to stderr)
   -Xmx*          :: a memory specification passed to java (default: -Xmx1g)

In the following, the period-delimited prefix of the input file name
is referred to as PREFIX.

If the input file is already log-transformed, this script will move it
to OUTDIR/PREFIX.final.pcl; otherwise, the script will log-transform
the file, producing OUTDIR/PREFIX.final.pcl

PREFIX is also used to find a match in the info file based on the
DatasetID field.

The pcl filename is normally prefixed with "GDSID." or "GDSID_SETID."
where GDSID is the GEO id, for example "GSE217.", and "SETID" is a
character string that does not contain a period.

Example:
  #{BN} -Xmx2g GSE13219.averaged.pcl GSE13219.info = .

Version: 2011.06.15
EOH
  exit 0
end

### FUNCTIONS
def verbose(*msgs)
  if @@verbosity
    msgs.each { |msg| $stderr.puts msg }
  end
end

#Important configuration columns
GDSID_col = 1
org_col   = 2  # (actually no longer relevant)
log_col  = 21

# Ensure the paths end with a /
JARDIR = (ARGV[2] == "=") ?
     %x(cd "#{File.dirname($0)}"; pwd).chomp! + "/" :
         ARGV[2].end_with?("/") ? ARGV[2] : (ARGV[2] + "/")

OUTDIR = ARGV[3].end_with?("/") ? ARGV[3] : (ARGV[3] + "/") if ARGV[3]


#Read the configuration file and locate the needed info
filename = File.basename(ARGV[0])

i = filename.index(".")
if i
  GDSID = filename.slice(0, i)
else
  $stderr.print "#{BN}: the pcl filename must have a period-delimited prefix.\n"
  exit 1
end

logged = org = nil

# WARNING: The combination of Ruby1.9 and UTF8 requires care!
# Ruby1.9: File.open(ARGV[1], :encoding => "UTF-8").each_line do |line|
# Ruby1.8: IO.foreach(ARGV[1]) do |line|
begin
  fd = File.open(ARGV[1], :encoding => "UTF-8")
rescue
  fd = File.open(ARGV[1] )
end

fd.each_line do |line|
	line.chomp!
	parts = line.split("\t")
	if parts[GDSID_col] == GDSID
		org = parts[org_col].gsub(" ","_")
		logged  = parts[log_col].to_i
	end
end

#Make sure the organism directory exists
#if !File.exists?(OUTDIR + org)
#	Dir.mkdir(OUTDIR + org)
#end

#If the data needs to be log transformed, do it and pipe to the final file
if logged == 0
	cmd = "java " + jopt + " -jar " + JARDIR + "DivLogNorm.jar " + ARGV[0] + " 0 1 0 >" + OUTDIR + GDSID + ".final.pcl"
#Otherwise, just move the file
else
	cmd = "mv " + ARGV[0] + " " + OUTDIR + GDSID + ".final.pcl"
end

puts cmd
system cmd
