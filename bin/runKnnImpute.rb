#!/usr/bin/env ruby
# Copyright (C) 2007 Matt Hibbs
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.


if ARGV.length < 2
  puts "ARGS:"
  puts "0 - pcl file to impute"
  puts "1 - configuration file (often allInfo.txt)"
  puts "2 - path to KNN impute program"
  puts "3 - output file"
  puts "This script will run KNNimpute on the pcl file.  However, KNNImputer"
  puts "does not perform correctly on raw single channel data, so these data"
  puts "are first log transformed, then imputed, then exponentiated." 
  exit 0
end


#Important configuration columns
GDSID_col = 1
numChan_col = 20
log_col = 21

#Read the configuration file and locate the needed info
if ARGV[0].rindex("/") != nil
  filename = ARGV[0].slice(ARGV[0].rindex("/") + 1, ARGV[0].length)
else
  filename = ARGV[0]
end
GDSID = filename.slice(0, filename.index("."))

numChan = logged = nil
IO.foreach(ARGV[1]) do |line|
	line.chomp!
	parts = line.split("\t")
	if parts[GDSID_col] == GDSID
		numChan = parts[numChan_col].to_i
		logged = parts[log_col].to_i
	end
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

#Run KNNImputer
cmd = ARGV[2] + "KNNImputer -i " + inFileName + " -o " + outFileName + " -k 10 -m 0.7 -d euclidean"
$stderr.puts cmd
system cmd

#If this file was log transformed, expoentiate it so that probe averaging isn't affected
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
	

