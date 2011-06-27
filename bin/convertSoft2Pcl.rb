#!/usr/bin/env ruby
# Copyright (C) 2007 Matt Hibbs
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.


if ARGV.length < 1
	$stderr.puts "ARGS:"
	$stderr.puts "0 - soft file to convert into pcl"
	$stderr.puts "1 - filename to output the pcl formatted file to"
	exit 0
end

def is_numeric?(x)
	true if Float(x) rescue false
end

tag_org = "!dataset_platform_organism = "
tag_plat = "!dataset_platform = "
tag_dset = "^DATASET = "
tag_title = "!dataset_title = "
tag_desc = "!dataset_description = "
tag_pmid = "!dataset_pubmed_id = "
tag_featCount = "!dataset_feature_count = "
tag_chanCount = "!dataset_channel_count = "
tag_sampCount = "!dataset_sample_count = "
tag_valType = "!dataset_value_type = "
tag_modDate = "!dataset_update_date = "

tag_tableStart = "!dataset_table_begin"
tag_tableEnd = "!dataset_table_end"

labels = Hash.new
in_data_section = false
at_labels = false

min = 1e10
max = -1e10
mean = numPos = numNeg = numZero = numMissing = numTotal = 0

fout = File.open(ARGV[1],"w")

print "File\tDatasetID\tOrganism\tPlatform\tValueType\t#channels\tTitle\tDesctiption\tPubMedID\t"
print "#features\t#samples\tdate\tMin\tMax\tMean\t#Neg\t#Pos\t#Zero\t#MV\t#Total\t#Channels\t"
puts  "logged\tzerosAreMVs\tMVcutoff"

org = platform = dset = title = desc = pmid = featCount = chanCount = sampCount = valType = date = ""
file = ARGV[0]
IO.foreach(file) do |line|
	line.strip!
	
	if line.include?(tag_org)
		org = line.slice(tag_org.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_plat)
		platform = line.slice(tag_plat.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_dset)
		dset = line.slice(tag_dset.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_title)
		title = line.slice(tag_title.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_desc)
		desc = line.slice(tag_desc.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_pmid)
		pmid = line.slice(tag_pmid.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_featCount)
		featCount = line.slice(tag_featCount.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_chanCount)
		chanCount = line.slice(tag_chanCount.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_sampCount)
		sampCount = line.slice(tag_sampCount.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_valType)
		valType = line.slice(tag_valType.length, line.length).gsub("\t"," ")
	elsif line.include?(tag_modDate)
		date = line.slice(tag_modDate.length, line.length).gsub("\t"," ")
	end
	
	
	if (!in_data_section && line.slice(0,1) == "#")
		parts = line.split(" = ")
		key = parts[0].slice(1,parts[0].length).gsub("\t"," ").gsub("~","-")
		labels[key] = parts[1]
	end
	
	if at_labels
		at_labels = false
		in_data_section = true
		parts = line.split("\t")
		fout.print labels[parts[0]] + "\t" + labels[parts[1]] + "\tGWEIGHT"
		for i in 2...parts.length
			fout.print "\t" + labels[parts[i]]
		end
		fout.puts
		fout.print "EWEIGHT\t\t"
		for i in 2...parts.length
			fout.print "\t1"
		end
		fout.puts
	elsif line.include?(tag_tableStart)
		at_labels = true
	elsif in_data_section
		if line.include?(tag_tableEnd)
			in_data_section = false
		else
			parts = line.split("\t")
			fout.print parts[0] + "\t" + parts[1] + "\t1"
			for i in 2...parts.length
				if !is_numeric?(parts[i])
					numMissing += 1
					fout.print "\t"
				else
					numTotal += 1
					fout.print "\t" + parts[i]
					val = Float(parts[i])
					mean += val
					if val > max
						max = val
					end
					if val < min
						min = val
					end
					if val > 0
						numPos += 1
					elsif val < 0
						numNeg += 1
					else
						numZero += 1
					end
				end
			end
			fout.puts
		end
	end
end

tested_numChans = 2
tested_logXformed = 1
tested_zerosMVs = 0
tested_MVcutoff = "NA"

# LA: avoid division by zero error
if numTotal != 0
  mean /= numTotal
else
  mean = 0
end


if numNeg == 0 || (numPos / numNeg.to_f) > 7.5
	tested_numChans = 1
end
if max > 100
	tested_logXformed = 0
end
if numZero > 5*numMissing
	tested_zerosMVs = 1
end
if tested_numChans == 1
	if tested_logXformed == 1
		tested_MVcutoff = 0
	elsif min > -500
		tested_MVcutoff = 2
	end
end


fout.flush
fout.close

print file + "\t" + dset + "\t" + org + "\t" + platform + "\t" + valType + "\t" + chanCount + "\t" + title
print "\t" + desc + "\t" + pmid + "\t" + featCount + "\t" + sampCount + "\t" + date

print "\t" + min.to_s + "\t" + max.to_s + "\t" + mean.to_s + "\t" + numNeg.to_s + "\t" + numPos.to_s + "\t"
print numZero.to_s + "\t" + numMissing.to_s + "\t" + numTotal.to_s + "\t" + tested_numChans.to_s + "\t"
puts  tested_logXformed.to_s + "\t" + tested_zerosMVs.to_s + "\t" + tested_MVcutoff.to_s



