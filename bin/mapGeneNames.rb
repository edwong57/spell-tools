#!/usr/bin/env ruby
# Copyright (C) 2007 Matt Hibbs and (C) 2011 Peter Koppstein
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.

# This is the Troilkatt version as of May 27, 2011 with changes by peak@princeton.edu including:
#  - modifications to allow use with UTF-8 input using either Ruby 1.8 or 1.9
#  - map=PATHNAME as alternative to "organism file"

require 'set'

BN = File.basename($0)

case ARGV[0]
when "-v", "--verbose"
  @@verbosity = true
  ARGV.shift
else
  @@verbosity = false
end

if ARGV.length < 3
  print <<-EOH
Syntax: #{BN} [-v | --verbose] ARGS
ARGS:
   0 - input pcl filename, the period-delimited prefix of which is used to find a match in the info file (see below)
   1 - info file (aka configuration file - often GSEnnn.info or allInfo.txt)
   2 - the pathname of the organism file, or "map=" followed immediately by the pathname of the map file (see below)
   3 - output file

Options:
   -v | --verbose :: verbose (to stderr)

This script creates a new pcl file by copying the input pcl file after mapping the gene names in it to standard symbols.

The period-delimited prefix of the pcl filename is used to find a match in the info file based on the DatasetID field.
The pcl filename is normally prefixed with "GDSID." or "GDSID_SETID." where GDSID is the GEO id, for example "GSE217.",
and "SETID" is a character string that does not contain a period.

If specified, the organism file should consist of rows of the form

ORGANISM<tab>PATHNAME

where:
   ORGANISM is the organism name (e.g. 'Saccharomyces cerevisiae') as it appears in the .info file for this GDSID
   PATHNAME is the pathname of the alias map file for ORGANISM

The alias map file should consist of rows of the form:
  G S1|S2|....|Sn

where G is a gene name and the S are the corresponding standard symbols.
(Each input row for G results in n output rows.)

The .info file is consulted for the value of ORGANISM for this GDSID.

Example:
  #{BN} GSE13219.knn.pcl GSE13219.sfp.info map=sce.map

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

### VARIABLES
#Important configuration column indices (IO=0):
GDSID_col = 1
org_col = 2

###Read the configuration file and locate the needed info
filename = File.basename(ARGV[0])

i = filename.index(".")
if i
  GDSID = filename.slice(0, i)
else
  $stderr.print "#{BN}: the pcl filename must have a period-delimited prefix.\n"
  exit 1
end

verbose("GDSID=#{GDSID}")

org = nil

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
        # verbose(line)
	parts = line.split("\t")
	if parts[GDSID_col] == GDSID
		org = parts[org_col]
	end
end

verbose("org=#{org}")

#Determine the mapping file:

orgs = Hash.new
if ARGV[2].start_with?("map=")
  orgs[org] = ARGV[2].sub("map=","")
  verbose("orgs[#{org}]=#{orgs[org]}")
else
  IO.foreach(ARGV[2]) do |line|
	line.chomp!
	if line.slice(0,1) != "#"
		parts = line.split("\t")
		orgs[parts[0]] = parts[1]
	end
  end
end

#If the organism is included in the list of mappings
if org != nil && orgs[org] != nil
	#Load the alias file
	map = Hash.new
	IO.foreach(orgs[org]) do |line|
		line.chomp!
		parts = line.upcase.split("\t")
		if parts.length > 1
			subparts = parts[1].split("|")
			if map[parts[0]] == nil
				map[parts[0]] = Set.new
			end
			map[parts[0]].merge(subparts)
		end
	end
	
	#Perform the mapping using the aliases provided
	fout = File.open(ARGV[3],"w")
	IO.foreach(ARGV[0]) do |line|
		line.chomp!
		if $. < 3
			fout.puts line
		else
			parts = line.upcase.split("\t")
			if map[parts[1]] != nil
				for name in map[parts[1]].to_a.sort
					fout.print name + "\t" + name
					for i in 2...parts.length
						fout.print "\t" + parts[i]
					end
					fout.puts
				end
			else
				$stderr.puts "WARNING: " + parts[1] + " could not be mapped"
			end
		end
	end
	fout.flush
	fout.close
end

