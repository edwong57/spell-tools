#!/usr/bin/ruby
# Copyright (C) 2007 Matt Hibbs
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.


require 'set'

if ARGV.length < 3
	puts "ARGS:"
	puts "0 - pcl file to map aliases from"
	puts "1 - config file (often allInfo.txt)"
	puts "2 - organism file, containing name in column 1, path to alias map in column 2"
  puts "3 - output file"
	puts "Creates a new file with \".map\" appended that contains a new pcl"
	puts "file with gene names mapped to standard symbols"
	exit 0
end

#Read the organism information file
orgs = Hash.new
IO.foreach(ARGV[2]) do |line|
	line.chomp!
	if line.slice(0,1) != "#"
		parts = line.split("\t")
		orgs[parts[0]] = parts[1]
	end
end

#Important configuration columns
GDSID_col = 1
org_col = 2

#Read the configuration file and locate the needed info
if ARGV[0].rindex("/") != nil
  filename = ARGV[0].slice(ARGV[0].rindex("/") + 1, ARGV[0].length)
else
  filename = ARGV[0]
end
GDSID = filename.slice(0, filename.index("."))

org = nil
IO.foreach(ARGV[1]) do |line|
	line.chomp!
	parts = line.split("\t")
	if parts[GDSID_col] == GDSID
		org = parts[org_col]
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

