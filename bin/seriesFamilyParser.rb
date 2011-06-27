#!/usr/bin/env ruby
# Copyright (C) 2007-2011 Matt Hibbs, Peter Koppstein
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include this copyright notice.

# Source: troilkatt/scripts/array_utils/seriesFamilyParser.rb as of Feb 4. 2011
# (Except for white space, the above version is identical to the version provided by Patrick Bradley (11 Jan 2011).)
# This version contains modifications by peak@princeton.edu marked #peak hereafter.

require 'set'

#peak:
require 'optparse'
require 'ostruct'

# Unless otherwise specified, for any key, options.key is nil
options = OpenStruct.new

options.version = "0.0.1"

optionparser = OptionParser.new do |opts|
  bn = File.basename($0)
  opts.banner = "Usage: #{bn} [options] FILE [configurationFile]"

  opts.separator("FILE should be an NCBI GSE \"soft\" file.")
  opts.separator("This script emits data in PCL format to stdout.")
  opts.separator("Options:")

  opts.on("-b", "--batch", "Batch mode (no interactions)") do
    options.batch = true
  end

  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    options.verbose = v
  end

  opts.on_tail("-V", "--version", "Show version and exit") do
    puts options.version
    exit
  end

  explanation = <<-EOS
configurationFile:
  If configurationFile is specified, the first four lines should specify the column numbers for:
    Platform ProbeID
    Platform GeneName
    Sample ProbeID
    Sample Value

    The next line should specify 1 for convertZerosToMV (convert zeros to missing values), else 0.

Example:
    #{bn} --batch GSE0100_family.soft > GSE0100.pcl

EOS

  opts.on("-h", "--help", "print this help") do
    puts opts
    print "\n", explanation
    exit 0
  end
end

# Modify ARGV
optionparser.parse!


#peak: message (a string) is optional:
def cancelout(message = "Received negative number; canceling...", rc=0)
  $stderr.print(message + "\n" ) if message
  Process.exit(rc)
end


if ARGV.length < 1
  #peak: 
  # To emit help: optionparser.parse("-h")
  puts optionparser.banner
  puts "For help type #{File.basename($0)} -h"
  exit 0
end

inPlatTable = false
inSampTable = false
inSample = false

platIdIdxs = Hash.new			#Per-platform, keep which column has unique IDs
platNameIdxs = Hash.new		#Per-platform, keep which column has gene names

sampIdIdxs = Hash.new 		#Per-platform, which column of the sample table has unique IDs
sampValueIdxs = Hash.new 	#Per-platform, which column of the sample table has desired value

sampPlatform = Hash.new		#per-sample, which platform to use
sampTitle = Hash.new			#per-sample, title string for sample

platLines = nil
sampLines = nil
platID = nil
sampID = nil

missingValueCount = 0
zeroValueCount = 0

class Array
  # return the index in self of the first element in elts that occurs in self; else nil
  def first_match(elts)
    hash = {}
    if self.size == 0
      return nil
    end
    for i in 0 ... self.size
      hash[self[i]] = i
    end 
    for i in 0 ... elts.size
      x = hash[elts[i]]
      if x != nil
        return x
      end
    end
    return nil
  end
end

# peak -end.

#Adding config script hack

if ARGV.length ==2
  $stderr.print("Will use parameters from config file\n")
  useConfig=true
  file = File.new(ARGV[1], "r")
  line=file.gets.to_i
  platIdIdxAll=line

  line=file.gets.to_i
  platNameIdxAll=line

  line=file.gets.to_i
  sampIdIdxAll=line

  line=file.gets.to_i
  sampValueIdxAll=line

  line=file.gets.to_i
  convertZerosToMV=line

  $stderr.printf("Platform ProbeId is in column %d\n", platIdIdxAll)
  $stderr.printf("Platform GeneName is in column %d\n", platNameIdxAll)
  $stderr.printf("Sample ProbeId is in column %d\n", sampIdIdxAll)
  $stderr.printf("Sample Value is in column %d\n", sampValueIdxAll)
  if convertZerosToMV==1
    $stderr.print("Will convert zeros to missing values\n")
  else
    $stderr.print("Will NOT convert zeros to missing values\n")
  end
else
  useConfig=false
end

#Do an initial pass through to determine information about platforms, samples, etc.
IO.foreach(ARGV[0]) do |line|
  line.chomp!

  if line.include?("^PLATFORM")
    platID = line.slice(12,line.length)
    platLines = Array.new
  end
  if line.include?("!platform_table_begin")
    inPlatTable = true
  elsif inPlatTable == true
    if line.include?("!platform_table_end")
      inPlatTable = false
    elsif platLines.length < 5
      platLines.push(line)
    end
    if platLines.length == 5 && platIdIdxs[platID] == nil
      if useConfig
        platIdIdxs[platID] = platIdIdxAll
        platNameIdxs[platID] = platNameIdxAll

      #peak - start
      elsif options.batch
        # Platform probe id
        headers = platLines[0].split("\t")
        j = headers.first_match([ "ID" ])
        if j == nil: cancelout("unable to find label for platform probe id", 1) end
        platIdIdxs[platID] = j
        $stderr.print("Platform probe id is being set to " + j.to_s + " (" + headers[j] + ")\n")

        # Gene names
        # "ORF", "GB_ACC", "GENOME_ACC", "RANGE_GB", "GB_LIST", "Gene_ID", "CLONE_ID", "Genename", "GENE_NAME", "Common Name", "Gene Symbol", "GENE_SYMBOL"
        # ordered list of candidates from SeriesFamilyParser.java:
        j = headers.first_match( [ "GENE_SYMBOL", "GENE SYMBOL", "GENE_NAME", "GENE NAME", "GENE ID", "ORF", "DDB" ] )
        if j == nil: cancelout("unable to find label for gene names", 1) end
        platNameIdxs[platID] = j
        $stderr.print("gene name column is being set to " + j.to_s + " (" + headers[j] + ")\n")
      #peak - end

      else
        $stderr.puts "Columns for platform " + platID
        numCols = platLines[0].split("\t").length
        for j in 0...numCols
          $stderr.print j.to_s
          for i in 0...platLines.length
            parts = platLines[i].split("\t")
            if parts.length > j
              $stderr.print "\t" + parts[j]
            else
              $stderr.print "\t--"
            end
          end
          $stderr.puts
        end
        $stderr.print "Which row contains the probe ID? "
        platIdIdxs[platID] = $stdin.gets.to_i
        # useful if called from batch file -- pjbradle
        # ADD END, USEFUL MSG
        if platIdIdxs[platID] < 0: cancelout() end
        $stderr.print "Which row contains the gene names? "
        platNameIdxs[platID] = $stdin.gets.to_i
        if platNameIdxs[platID] < 0: cancelout() end
      end
    end
  end

  if line.include?("^SAMPLE")
    inSample = true
    sampID = line.slice(10,line.length)
  elsif inSample == true
    if line.include?("!Sample_title")
      sampTitle[sampID] = line.slice(16,line.length)
    elsif line.include?("!Sample_platform_id")
      sampPlatform[sampID] = line.slice(22,line.length)
    elsif line.include?("!sample_table_begin")
      inSampTable = true
      sampLines = Array.new
    elsif inSampTable == true
      if line.include?("!sample_table_end")
        inSampTable = false
        inSample = false
      elsif sampLines.length < 5
        sampLines.push(line)
      end
      if sampLines.length == 5 && sampIdIdxs[sampPlatform[sampID]] == nil
        if useConfig
          sampIdIdxs[sampPlatform[sampID]] = sampIdIdxAll
          sampValueIdxs[sampPlatform[sampID]] = sampValueIdxAll
        #peak - start
        elsif options.batch
          headers=sampLines[0].split("\t")
          #!! $stderr.puts "Headers", headers, "\n"

          # List of keywords (most-preferred first) used to detect the sample probe ID
          j = headers.first_match( ["ID_REF"] )
          if j == nil: cancelout("field for sample probeID not found",1 ) end
          sampIdIdxs[sampPlatform[sampID]] = j
          $stderr.print("Sample probe ID column is being set to " + j.to_s + " (" + headers[j] + ")\n")

          # List of keywords (most-preferred first) used to detect the expression values:
          j = headers.first_match( ["VALUE"] )
          if j == nil: cancelout("field for expression value not found",1 ) end
          sampValueIdxs[sampPlatform[sampID]] = j
          $stderr.print("Sample probe ID column is being set to " + j.to_s + " (" + headers[j] + ")\n")
        #peak - end
        else
          $stderr.puts "Columns for sample " + sampID + " using platform " + sampPlatform[sampID].to_s
          numCols = sampLines[0].split("\t").length
          for j in 0...numCols
            $stderr.print j.to_s
            for i in 0...sampLines.length
              parts = sampLines[i].split("\t")
              if parts.length > j
                $stderr.print "\t" + parts[j]
              else
                $stderr.print "\t--"
              end
            end
            $stderr.puts
          end

          $stderr.print "Which row contains the probe ID? "
          sampIdIdxs[sampPlatform[sampID]] = $stdin.gets.to_i
          if sampIdIdxs[sampPlatform[sampID]] < 0: cancelout() end
          $stderr.print "Which row contains the expression values? "
          sampValueIdxs[sampPlatform[sampID]] = $stdin.gets.to_i
          if sampValueIdxs[sampPlatform[sampID]] < 0: cancelout() end
        end
      end
      if sampIdIdxs[sampPlatform[sampID]] != nil
        parts = line.split("\t")
        val = parts[sampValueIdxs[sampPlatform[sampID]]]
        if val == ""
          missingValueCount += 1
        elsif val.to_f == 0
          zeroValueCount += 1
        end
      end
    end
  end
end
if useConfig
  zerosAsMissingVals=convertZerosToMV
elsif options.batch && missingValueCount > 0
  $stderr.puts "missingValueCount == " + missingValueCount.to_s + " > 0 so we will NOT treat exact zeros as missing values\n"
  zerosAsMissingVals = 0
else
  $stderr.puts "\nThere are " + missingValueCount.to_s + " missing values, and " + zeroValueCount.to_s + " exactly zero values"
  $stderr.print "Treat exact zeros as missing values (0=no, 1=yes)? "
  zerosAsMissingVals = $stdin.gets.to_i
end
$stderr.puts
$stderr.puts "Parsing file into pcl format..."

samples = Hash.new	#Hash from sampID to a Hash from geneName to value
platMaps = Hash.new	#Hash from platID to a Hash from uniqueID to geneName
idToName = nil
plat = nil
genes = Set.new
first = false

#On a second pass through, actually parse out everything
IO.foreach(ARGV[0]) do |line|
  line.chomp!

  if line.include?("^PLATFORM")
    platID = line.slice(12,line.length)
    platMaps[platID] = Hash.new
  end
  if line.include?("!platform_table_begin")
    inPlatTable = true
  elsif inPlatTable == true
    if line.include?("!platform_table_end")
      inPlatTable = false
    else
      parts = line.split("\t")
      platMaps[platID][parts[platIdIdxs[platID]]] = parts[platNameIdxs[platID]]
    end
  end

  if line.include?("^SAMPLE")
    inSample = true
    sampID = line.slice(10,line.length)
    samples[sampID] = Hash.new
    plat = sampPlatform[sampID]
    idToName = platMaps[plat]
  elsif inSample == true
    if line.include?("!sample_table_begin")
      inSampTable = true
      first = true
    elsif inSampTable == true
      if first == true
        first = false
      elsif line.include?("!sample_table_end")
        inSampTable = false
        inSample = false
        plat = nil
        idToName = nil
      else
        parts = line.split("\t")
        gname = idToName[parts[sampIdIdxs[plat]]]
        if gname != nil && gname != ""
          samples[sampID][gname] = parts[sampValueIdxs[plat]]
          genes.add(gname)
        else
          #$stderr.puts "WARNING: For sample " + sampID + ", " + parts[sampIdIdxs[plat]] + " was not present in the " + sampPlatform[sampID] + " platform table"
        end
      end
    end
  end
end

#Output the pcl file
sampArr = samples.keys.to_a.sort
print "YORF\tNAME\tGWEIGHT"
for i in 0...sampArr.length
  print "\t" + sampTitle[sampArr[i]]
end
puts
print "EWEIGHT\t\t"
for i in 0...sampArr.length
  print "\t1"
end
puts

for gene in genes.to_a.sort
  print gene + "\t" + gene + "\t1"
  for i in 0...sampArr.length
    val = samples[sampArr[i]][gene]
    if val == nil
      print "\t"
    else
      if zerosAsMissingVals == 1 && val.to_f == 0
        print "\t"
      else
        print "\t" + val
      end
    end
  end
  puts
end
