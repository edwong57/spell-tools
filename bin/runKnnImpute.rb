#!/usr/bin/env ruby
# Copyright (C) 2007 Matt Hibbs
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.

# This is the Troilkatt version as of May 27, 2011 with minor tweaks by peak@princeton.edu 

BN = File.basename($0)

case ARGV[0]
when "-v", "--verbose"
  verbose = true
  ARGV.shift
else
  verbose = false
end

if ARGV.length < 2
  print <<-EOH
Syntax: #{BN} [-v | --verbose] ARGS
ARGS:
   0 - pcl file to impute - pathname or filename
   1 - info file (aka configuration file - often GSEnnn.info or allInfo.txt)
   2 - path to KNN impute program
   3 - output filename

Options:
   -v | --verbose :: verbose (to stderr)

This script will run KNNimpute on the pcl file.  Since KNNImputer
does not perform correctly on raw single channel data, such data
are first log transformed, then imputed, then exponentiated. 

If the configuration file is not in "standard format" then the filename should
have the GDSID as its period-delimited prefix.

Example:
  #{BN} GSE13219.sfp.pcl GSE13219.sfp.info /usr/local/bin GSE13219.knn.pcl

Version: 2011.06.02
EOH
  exit 0
end


### FUNCTIONS
def linecount(filename)
  # count = %x{wc -l < #{filename}}.to_i
  count = File.foreach(filename).count
end

### VARIABLES
#Important configuration column indices (IO=0):
GDSID_col = 1
numChan_col = 20
log_col = 21

### Read the configuration file and locate the needed info

# First try: assume the info file is in "standard format"

filename = File.basename(ARGV[0])
GDSID = found = standard = numChan = logged = nil

# If GDSID cannot be discerned from the filename, then look in the info file directly if it has just two lines:

i=filename.index("_") || filename.index(".")

if i
  GDSID = filename.slice(0, i)
elsif 2 == linecount(ARGV[1])
  row=0
  IO.foreach(ARGV[1]) do |line|
    row += 1
    line.chomp!
    parts = line.split("\t")
    if row == 1 and parts[GDSID_col] == "DatasetID"
      standard = 1
      if parts[numChan_col] == "#Channels" and parts[log_col] == "logged"
        standard = 2
      elsif verbose
        $stderr.print "parts[GDSID_col]:", parts[GDSID_col], "\n"
        $stderr.print "parts[numChan_col]:", parts[numChan_col], "\n"
        $stderr.print "parts[log_col]:", parts[log_col], "\n"
      end
    end
    if row == 2 and standard
      GDSID=parts[GDSID_col]
      numChan = parts[numChan_col].to_i
      logged = parts[log_col].to_i
      found = 1
    end
  end
end

unless GDSID
  $stderr.puts "#{BN}: unable to determine GDSID from filename or info file"
  exit 1
end

if !found
  if verbose
    $stderr.puts "Searching for GDSID=#{GDSID} based on filename"
  end
  IO.foreach(ARGV[1]) do |line|
    line.chomp!
    parts = line.split("\t")
    if parts[GDSID_col] == GDSID
      found = 2
      numChan = parts[numChan_col].to_i
      logged = parts[log_col].to_i
      break
    end
  end
end

if !found or !numChan or !logged
  $stderr.puts "#{BN}: unable to make sense of info file #{ARGV[1]}:"
  case found
  when 1
    $stderr.puts "GDSID in info file is #{GDSID}."
  when 2
    $stderr.puts "GDSID, which was computed from filename as #{GDSID}, was found in info file."
  else
    $stderr.puts "GDSID, which was computed from filename as #{GDSID}, was not found in info file."
  end
  exit 1
end
  

inFileName = ARGV[0]
outFileName = ARGV[3]

#If this file needs to be log transformed, do it prior to imputing
if logged == 0
	#Create a .tmp file that contains the log xformed version
	fout = File.open(ARGV[3] + ".tmp1","w")
	IO.foreach(ARGV[0]) do |line|
		line.chomp!
		if $. < 3
			fout.puts line
		else
			parts = line.split("\t")
			fout.print parts[0] + "\t" + parts[1] + "\t" + parts[2]
			for i in 3...parts.length
				begin
					val = Float(parts[i])
					fout.print "\t" + Math.log(val).to_s
				rescue
					fout.print "\t"
				end
			end
			fout.puts
		end
	end
	fout.flush
	fout.close
	#Set the in and out file names appropriately
	inFileName = ARGV[3] + ".tmp1"
	outFileName = ARGV[3] + ".tmp2"
end

# Run KNNImputer
# Ensure we have just one /
cmd = ARGV[2]
cmd="${cmd}/" unless cmd[-1] == "/"[0]

cmd += "KNNImputer -i " + inFileName + " -o " + outFileName + " -k 10 -m 0.7 -d euclidean"
$stderr.puts cmd
system cmd

#If this file was log transformed, exponentiate it so that probe averaging isn't affected
if logged == 0
	fout = File.open(ARGV[3], "w")
	IO.foreach(outFileName) do |line|
		line.chomp!
		if $. < 3
			fout.puts line
		else
			parts = line.split("\t")
			fout.print parts[0] + "\t" + parts[1] + "\t" + parts[2]
			for i in 3...parts.length
				fout.print "\t" + Math.exp(parts[i].to_f).to_s
			end
			fout.puts
		end
	end
	fout.flush
	fout.close
	#Delete the intermediate files
	cmd = "rm " + inFileName
	#$stderr.puts cmd
	system cmd
	cmd = "rm " + outFileName
	#$stderr.puts cmd
	system cmd
end
