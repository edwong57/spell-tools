#!/usr/bin/ruby
# Copyright (C) 2007 Matt Hibbs
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.


if ARGV.length < 3
	puts "ARGS:"
	puts "0 - pcl file"
	puts "1 - configuration file (often allInfo.txt)"
	puts "2 - path to DivLogNorm.jar"
	puts "3 - path to final data folders"
	puts "Script will log transform the file if needed and move it to an"
	puts "organism specific sub-directory"
	exit 0
end

#Important configuration columns
GDSID_col = 1
org_col = 2
log_col = 21

#Read the configuration file and locate the needed info
if ARGV[0].rindex("/") != nil
  filename = ARGV[0].slice(ARGV[0].rindex("/") + 1, ARGV[0].length)
else
  filename = ARGV[0]
end
GDSID = filename.slice(0, filename.index("."))

logged = org = nil
IO.foreach(ARGV[1]) do |line|
	line.chomp!
	parts = line.split("\t")
	if parts[GDSID_col] == GDSID
		org = parts[org_col].gsub(" ","_")
		logged  = parts[log_col].to_i
	end
end

#Make sure the organism directory exists
#if !File.exists?(ARGV[3] + org)
#	Dir.mkdir(ARGV[3] + org)
#end

#If the data needs to be log transformed, do it and pipe to the final file
if logged == 0
	cmd = "java -Xmx1g -jar " + ARGV[2] + "DivLogNorm.jar " + ARGV[0] + " 0 1 0 >" + ARGV[3] + GDSID + ".final.pcl"
	puts cmd
	system cmd
#Otherwise, just move the file
else
	cmd = "mv " + ARGV[0] + " " + ARGV[3] + GDSID + ".final.pcl"
	puts cmd
	system cmd
end


